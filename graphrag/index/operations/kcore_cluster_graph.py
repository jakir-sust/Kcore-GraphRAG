from pickle import TRUE
import networkx as nx
import logging
from itertools import combinations

from collections import deque, Counter, defaultdict
from typing import List, Tuple
from graphrag.index.utils.stable_lcc import stable_largest_connected_component
Communities = List[Tuple[int, int, int, list[str]]]

logger = logging.getLogger(__name__)

from graspologic.partition import leiden
from collections import defaultdict, deque
import networkx as nx
import matplotlib.pyplot as plt
import networkx as nx
import plotly.graph_objects as go

def split_connected_chunks(G, nodes, max_size=20):
    nodes = set(nodes)
    tmp_nodes = nodes.copy()
    chunks = []

    subG = G.subgraph(nodes)
    deg = dict(subG.degree())

    singleton_list= []
    node2chunk = {}


    while nodes:
        # 1. pick densest seed
        seed = max(nodes, key=lambda x: (deg.get(x, 0), hash(x)))
        S = {seed}
        frontier = set(subG.neighbors(seed)) & nodes

        # track internal edges
        internal_edges = 0

        # 2. grow the set greedily
        while len(S) < max_size and frontier:
            best_node = None
            best_gain = -1

            for v in sorted(frontier):
                # number of edges v would add inside S
                gain = sum(1 for u in S if subG.has_edge(u, v))
                if gain > best_gain:
                    best_gain = gain
                    best_node = v

            if best_node is None:
                break

            S.add(best_node)
            frontier |= set(subG.neighbors(best_node)) & nodes
            frontier -= S

        chunk_id = len(chunks)
        for n in S:
            node2chunk[n] = chunk_id

        if len(S) == 1:
            singleton_list.extend(list(S))
        else:
            # print(len(list(S)))
            chunks.append(list(S)) 
        nodes -= S
            

    for s in singleton_list:
        neighbor_chunks = {node2chunk[v] for v in subG.neighbors(s) if v in node2chunk}
        if neighbor_chunks:
            # assign to first neighbor chunk (could also merge multiple if needed)
            target_chunk_id = list(neighbor_chunks)[0]
            chunks[target_chunk_id].append(s)   # <-- update the chunks list
            node2chunk[s] = target_chunk_id
    return chunks


def split_2hop_connected_chunks(G, nodes, anchor_nodes, max_size=20):
    nodes = set(nodes)
    anchor_nodes = set(anchor_nodes)
    chunks = []
    singleton_nodes = []        # NEW
    node2chunk = {}             # NEW
    chunk_id = 0                # NEW

    subG = G.subgraph(nodes | anchor_nodes)
    deg = dict(subG.degree())


    # Precompute anchor sets per node
    A = {u: set(G.neighbors(u)) & anchor_nodes for u in nodes}

    while nodes:
        # 1) pick seed with max anchor coverage
        seed = max(nodes, key=lambda u: (len(A[u]), hash(u)))
        S = {seed}

        # frontier: nodes sharing at least one anchor with seed
        frontier = {v for v in nodes if v != seed and A[v] & A[seed]}

        # 2) greedy growth
        while len(S) < max_size and frontier:
            best_node, best_gain = None, -1

            for v in frontier:
                gain = sum(len(A[v] & A[u]) for u in S)
                if gain > best_gain:
                    best_gain = gain
                    best_node = v

            if best_node is None or best_gain == 0:
                break

            S.add(best_node)
            frontier |= {v for v in nodes if A[v] & A[best_node]}
            frontier -= S

        # 3) minimal anchors: connect ≥2 nodes in S
        anchors = {
            a for a in anchor_nodes
            if sum(1 for u in S if G.has_edge(u, a)) >= 2
        }

        chunk_nodes = list(S) + list(anchors)
        chunks.append((chunk_id, chunk_nodes))
        for u in S:
            node2chunk[u] = chunk_id
        chunk_id += 1

        nodes -= S

    chunks = [chunk_nodes for _, chunk_nodes in chunks]
    return chunks


def assign_small_singletons_to_clusters(G, clusters, singletons, node2cluster):
    cluster_id_counter = max(cid for _, cid, _, _ in clusters)
    updated_clusters = clusters.copy()

    unassigned = {tuple(c) for c in singletons}
    left_out_nodes = []

    def cluster_neighbors(cluster):
        nbrs = set()
        for u in cluster:
            nbrs.update(G.neighbors(u))
        return nbrs

    def compute_score(cluster):
        return sum(
            1 for nbr in cluster_neighbors(cluster)
            if nbr in node2cluster
        )

    while unassigned:
        small_cluster = max(unassigned, key=compute_score)

        # ---- Count neighbors per cluster ----
        neighbor_counts = defaultdict(int)
        for neighbor in cluster_neighbors(small_cluster):
            if neighbor in node2cluster:
                neighbor_counts[node2cluster[neighbor]] += 1

        if neighbor_counts:
            best_cluster_id = max(
                neighbor_counts.items(), key=lambda x: x[1]
            )[0]

            for idx, (lvl, cid, pid, nodes) in enumerate(updated_clusters):
                if cid == best_cluster_id:
                    updated_clusters[idx] = (
                        lvl, cid, pid, nodes + list(small_cluster)  
                    )
                    break

            for u in small_cluster:
                node2cluster[u] = best_cluster_id

        else:
            print("New cluster: ", small_cluster)
            left_out_nodes.extend(small_cluster)

            cluster_id_counter += 1
            updated_clusters.append(
                (0, cluster_id_counter, -1, list(small_cluster))  
            )

            for u in small_cluster:
                node2cluster[u] = cluster_id_counter

        unassigned.remove(small_cluster)

    return updated_clusters


def kcore_cluster_graph(graph: nx.Graph, use_lcc: bool = True, cluster_type="RkH") -> Communities:
    print("kcore_cluster_graph called: ", cluster_type)
    logger.info("Cluster type: %s", cluster_type)
    if use_lcc:
        graph = stable_largest_connected_component(graph)
        print(f"#Nodes and #Edges in LCC: {graph.number_of_nodes()}, {graph.number_of_edges()}")
    
    if len(graph) == 0:
        logger.warning("Graph has no nodes")
        return []

    # Compute k-core number for each node
    copy_graph = nx.Graph()
    copy_graph.add_nodes_from(graph.nodes(data=True))
    copy_graph.add_edges_from(graph.edges(data=True))
    copy_graph.remove_edges_from(list(nx.selfloop_edges(copy_graph)))

    core_numbers = nx.core_number(copy_graph)

    # Group nodes by k-core number (level)
    levels = sorted(set(core_numbers.values()), reverse=False)  # higher k-core = higher level
    clusters: Communities = []

    cnt = dict(sorted(Counter(core_numbers.values()).items(), key = lambda x:x[0], reverse=False))
    print("levels: ", set(levels), cnt)


    cluster_queue = deque()
    cluster_queue.append((1, -1, list(copy_graph.nodes()), 0))

    MAX_CLUSTER_SIZE, MIN_CLUSTER_SIZE = 20, 1
    
    if cluster_type in {"M2hC", "MRC"}:
        MIN_CLUSTER_SIZE = 3

    level = 1
    cluster_id_counter = 1

    # Keep track of node -> cluster_id and attach singleton
    node2cluster = {}  #
    global_singletons, global_small_clusters = [], []

    total_small_cluster = 0
   
    while level <= max(levels): 
        print("\n", "***" * 10, "Level= ", level, "***" * 10)
        cur_level_clusters = []
        q_size = len(cluster_queue)
        singletons, singletons_nodes = 0, []

        tmp_cluster_counter = cluster_id_counter
        two_hop_sizes = []
        for i in range(q_size):
            cluster_id, parent_id, cur_nodes, is_res = cluster_queue.popleft()

            nodes_in_level = [n for n in cur_nodes if core_numbers[n] >= level]
            cur_res_nodes = [n for n in cur_nodes if core_numbers[n] < level]
            
            print(f"level: {level}, q-index = {i},  Level_node = {len(nodes_in_level)},  Res_node= {len(cur_res_nodes)}")

            # Create a subgraph for this k-core level nodes
            nodes_in_level_clusters = []
            subgraph = copy_graph.subgraph(nodes_in_level)
            for component in nx.connected_components(subgraph):
                child_nodes = list(component)
                if level < max(levels) or len(child_nodes) <= MAX_CLUSTER_SIZE:
                    cluster_id_counter += 1; nodes_in_level_clusters.append(len(child_nodes))
                    clusters.append((level, cluster_id_counter, cluster_id, child_nodes))
                    node2cluster.update({n: cluster_id_counter for n in child_nodes})
                    cluster_queue.append((cluster_id_counter, cluster_id, child_nodes, 0))
                elif level == max(levels):
                    subclusters = split_connected_chunks(copy_graph, child_nodes, max_size=MAX_CLUSTER_SIZE)
                    for sub_nodes in subclusters:
                        # PRINT_SINGLETON += 1 if len(sub_nodes) == 1 else 0
                        cluster_id_counter += 1; nodes_in_level_clusters.append(len(sub_nodes))
                        clusters.append((level, cluster_id_counter, cluster_id, sub_nodes))
                        node2cluster.update({n: cluster_id_counter for n in sub_nodes})
            
            # Subgraphs of residual nodes (k-1 shell nodes)
            res_nodes_clusters = []
            shell_subgraph = copy_graph.subgraph(cur_res_nodes)
            # Find connected components within the level-1
            for component in nx.connected_components(shell_subgraph):
                nodes = list(component)
                if len(nodes) == 1: 
                    singletons += 1; singletons_nodes.extend(nodes)
                elif cluster_type in {"RkH"}:
                    if len(nodes) > MAX_CLUSTER_SIZE:
                        subclusters = split_connected_chunks(copy_graph, nodes, max_size=MAX_CLUSTER_SIZE)
                        for sub_nodes in subclusters:
                            if len(sub_nodes) == 1:
                                singletons += 1; singletons_nodes.extend(sub_nodes)
                                continue
                            cluster_id_counter += 1; res_nodes_clusters.append(len(sub_nodes))
                            clusters.append((level, cluster_id_counter, cluster_id, sub_nodes))
                            node2cluster.update({n: cluster_id_counter for n in sub_nodes})
                    else: 
                        cluster_id_counter += 1; res_nodes_clusters.append(len(nodes))
                        clusters.append((level, cluster_id_counter, cluster_id, nodes))  
                        node2cluster.update({n: cluster_id_counter for n in nodes})
                elif cluster_type in {"M2hC", "MRC"}:
                    if cluster_type == "MRC":
                        MIN_RES_CC_SIZE = MIN_CLUSTER_SIZE
                    else:
                        MIN_RES_CC_SIZE = 2

                    if len(nodes) <= MAX_CLUSTER_SIZE and len(nodes) >= MIN_RES_CC_SIZE:
                        cluster_id_counter += 1; res_nodes_clusters.append(len(nodes))
                        clusters.append((level, cluster_id_counter, cluster_id, nodes))  
                        node2cluster.update({n: cluster_id_counter for n in nodes})
                    elif len(nodes) > MAX_CLUSTER_SIZE:
                        PRINT_SINGLETON = 0
                        subclusters = split_connected_chunks(copy_graph, nodes, max_size=MAX_CLUSTER_SIZE)
                        for sub_nodes in subclusters:
                            PRINT_SINGLETON += 1 if len(sub_nodes) == 1 else 0
                            if len(sub_nodes) < MIN_RES_CC_SIZE:
                                singletons += 1; singletons_nodes.extend(sub_nodes)
                                continue
                            cluster_id_counter += 1; res_nodes_clusters.append(len(sub_nodes))
                            clusters.append((level, cluster_id_counter, cluster_id, sub_nodes))
                            node2cluster.update({n: cluster_id_counter for n in sub_nodes})
                    elif len(nodes) < MIN_RES_CC_SIZE:
                         singletons += 1; singletons_nodes.extend(nodes)
                        
            # within shell: 2-hop CC
            res_2hop_CC_clusters = []
            actual_singleton = []
            visited = set()
            for i, singleton_node in enumerate(singletons_nodes):
                if singleton_node in visited:
                    continue
                visited.add(singleton_node)
                cur_subgraphs = [singleton_node]

                medium_nodes = set()  

                queue = deque([singleton_node])
                while queue:
                    u = queue.popleft()
                    for v in copy_graph[u]: # 1-hop
                        for w in copy_graph[v]:
                            if w != u and w not in visited and w in singletons_nodes: # 2-hop
                                visited.add(w)
                                cur_subgraphs.append(w)
                                medium_nodes.add(v)  
                if len(cur_subgraphs) == 1:
                    actual_singleton.append(cur_subgraphs)
                elif cluster_type in {"RkH"}:
                    if len(cur_subgraphs) <= MAX_CLUSTER_SIZE:
                        cur_subgraphs.extend(sorted(medium_nodes))
                        cluster_id_counter += 1; res_2hop_CC_clusters.append(len(cur_subgraphs))
                        clusters.append((level, cluster_id_counter, cluster_id, cur_subgraphs))  
                        node2cluster.update({n: cluster_id_counter for n in cur_subgraphs})
                    else:
                        subclusters = split_2hop_connected_chunks(copy_graph, cur_subgraphs, medium_nodes, max_size=MAX_CLUSTER_SIZE)
                        for sub_nodes in subclusters:
                            PRINT_SINGLETON += 1 if len(sub_nodes) == 1 else 0
                            if len(sub_nodes) == 1:
                                actual_singleton.append(sub_nodes)
                                continue
                            cluster_id_counter += 1;  res_2hop_CC_clusters.append(len(sub_nodes))
                            clusters.append((level, cluster_id_counter, cluster_id, sub_nodes))
                            node2cluster.update({n: cluster_id_counter for n in sub_nodes})

                elif cluster_type in{"M2hC", "MRC"}:
                    if len(cur_subgraphs) <= MAX_CLUSTER_SIZE and len(cur_subgraphs) >= MIN_CLUSTER_SIZE:
                        cur_subgraphs.extend(sorted(medium_nodes))
                        cluster_id_counter += 1; res_2hop_CC_clusters.append(len(cur_subgraphs))
                        clusters.append((level, cluster_id_counter, cluster_id, cur_subgraphs))  
                        node2cluster.update({n: cluster_id_counter for n in cur_subgraphs})
                    elif len(cur_subgraphs) > MAX_CLUSTER_SIZE:
                        subclusters = split_2hop_connected_chunks(copy_graph, cur_subgraphs, medium_nodes, max_size=MAX_CLUSTER_SIZE)
                        for sub_nodes in subclusters:
                            PRINT_SINGLETON += 1 if len(sub_nodes) == 1 else 0
                            if len(sub_nodes) < MIN_CLUSTER_SIZE:
                                actual_singleton.append(sub_nodes)
                                continue
                            cluster_id_counter += 1;  res_2hop_CC_clusters.append(len(sub_nodes))
                            clusters.append((level, cluster_id_counter, cluster_id, sub_nodes))
                            node2cluster.update({n: cluster_id_counter for n in sub_nodes})
                    elif len(cur_subgraphs) < MIN_CLUSTER_SIZE:
                         singletons += 1; actual_singleton.append(cur_subgraphs)
                        
            sm_size = 2
            total_small_cluster += sum(1 for cl_sz in nodes_in_level_clusters if cl_sz == sm_size) 
            total_small_cluster += sum(1 for cl_sz in res_nodes_clusters if cl_sz == sm_size) 
            total_small_cluster += sum(1 for cl_sz in res_2hop_CC_clusters if cl_sz == sm_size+1) 

            # after S2CC and 2-hop CC
            after_2hop_CC_clusters = []
            if len(actual_singleton) > 0:
                # print(f"SINGLETON AFTER TWO HOP: {len(actual_singleton)}")
                all_singleton_nodes = [v for u in actual_singleton for v in u]
                global_singletons.extend(all_singleton_nodes)
                for small_cluster in actual_singleton:
                    global_small_clusters.append(small_cluster)
        
        cur_cluster_cnt = cluster_id_counter - tmp_cluster_counter
        # print(f"#Child: {q_size}, level: {level}, #cluster: {cur_cluster_cnt}, #tmp_singleton: {singletons}, #final_singleton: {len(actual_singleton)}")
        level += 1
        nodes_in_level_clusters.sort(reverse=True); res_nodes_clusters.sort(reverse=True); 
        res_2hop_CC_clusters.sort(reverse=True); after_2hop_CC_clusters.sort(reverse=True); 
        
        cur_level_clusters = nodes_in_level_clusters + res_nodes_clusters + res_2hop_CC_clusters + after_2hop_CC_clusters 
        
        # if True:
        #     print(f"nodes_in_level= {nodes_in_level_clusters}, res_nodes_clusters = {res_nodes_clusters}")
        #     print(f"res_2hop_CC_clusters = {res_2hop_CC_clusters}")
        #     print(f"after_2hop_CC_clusters= {after_2hop_CC_clusters}")
        #     small_cluster_len = sorted([len(cl) for cl in global_small_clusters])
        #     print("global_singletons len: ", len(global_singletons), small_cluster_len)
        #     cur_level_clusters = sorted(cur_level_clusters, reverse=True)
        #     # print(f"LEVEL COMPLETED--> Size = {len(cur_level_clusters)},  All_clusters = {cur_level_clusters}")
        #     print(f"LEVEL {level-1} COMPLETED <--------> Size = {len(cur_level_clusters)}")

        #     logger.info("global_singletons len: %d", len(global_singletons))
        #     logger.info(f"LEVEL {level-1} COMPLETED <--------> Size = {len(cur_level_clusters)}")

    clusters = assign_small_singletons_to_clusters(copy_graph, clusters, global_small_clusters, node2cluster)

    print("Total Cluster found: ", len(clusters), total_small_cluster, total_small_cluster /len(clusters) * 100)
    print("####" * 10, "HIERARCHY COMPLETED", "####" * 10, "\n\n")

    level_counts = Counter(level for level, _, _, _ in clusters)
    level_counts_sorted = dict(sorted(level_counts.items(), key=lambda x: x[0]))
    # print(level_counts_sorted)

    logger.info("Total Clusters: %d", len(clusters))
    logger.info("Level Cluster count: %s", level_counts_sorted)

    return clusters

