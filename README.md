# Core-based Hierarchies for Efficient GraphRAG

This repository contains the code accompanying the paper  
**“Core-based Hierarchies for Efficient GraphRAG.”**

Our implementation builds directly on the official GraphRAG benchmarking framework released by Microsoft. We introduce k-core–based hierarchical community construction algorithms as drop-in replacements for Leiden-based community detection. The overall pipeline, datasets, and evaluation procedures remain unchanged.

## Code Base and Attribution

The code follows the same setup and execution procedure as the original GraphRAG benchmarking repository:

https://github.com/microsoft/graphrag-benchmarking-datasets

We retain:
- The same dataset format  
- The same indexing and query-time execution flow  
- The same evaluation and head-to-head comparison framework  

Our contributions are limited to modifications in the community construction and hierarchy generation components.

## Running the Code

The code is executed using the standard GraphRAG commands. No additional steps are required beyond those described in the original GraphRAG repository.

### Installing Necessary Packages

Run the following from the root directory:

```bash
pip install -e ./graphrag
````

### Prepare Your Workspace

Create a new directory to hold your input files (e.g., for the Kevin Scott podcast benchmark):

```bash
mkdir -p ./kevin_scott_podcasts/input
```

Place the unzipped input files into the `input` folder.

### Initialize the GraphRAG Workspace

```bash
python graphrag/cli/main.py init  --root kevin_scott_podcasts/
```

This command creates two files inside `./kevin_scott_podcasts/`:

* `.env`: contains the `GRAPHRAG_API_KEY` environment variable. Set it to your OpenAI or Azure OpenAI API key.
* `settings.yaml`: configures the GraphRAG pipeline. You may edit this file to customize pipeline behavior.

### Indexing

Indexing constructs the knowledge graph and hierarchical communities:

```bash
python graphrag/cli/main.py index --root ./kevin_scott_podcasts --community RkH
```

Here, `RkH` specifies the community construction algorithm. The default is the original GraphRAG Leiden algorithm. You may alternatively select one of our proposed heuristics:

* `RkH`
* `M2hC`
* `MRC`

Indexing time varies depending on dataset size and may take several minutes to several hours.

### Query

You can now run queries against your indexed dataset. For example, to perform global search on leaf-level (`LF`) communities:

```bash
python graphrag/cli/main.py query \
  --root ./kevin_scott_podcasts \
  --method global \
  --community-level LF \
  --query "What recurring topics do tech leaders emphasize in their discussions?"
```
