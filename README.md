# ThunderRW: An In-Memory Graph Random Walk Engine

## Introduction

As random walk is a powerful tool in many graph processing, mining and
learning applications, this paper proposes an efficient in-memory random
walk engine named ThunderRW. Compared with existing parallel systems on
improving the performance of a single graph operation, ThunderRW supports
massive parallel random walks. The core design of ThunderRW is motivated
by our profiling results: common RW algorithms have as high as 73.1% CPU
pipeline slots stalled due to irregular memory access, which suffers
significantly more memory stalls than the conventional graph workloads
such as BFS and SSSP. To improve the memory efficiency, we first design a
generic step-centric programming model named Gather-Move-Update to abstract
different RW algorithms. Based on the programming model, we develop the step
interleaving technique to hide memory access latency by switching the executions
of different random walk queries. In our experiments, we use four representative
RW algorithms including PPR, DeepWalk, Node2Vec and MetaPath to demonstrate the
efficiency and programming flexibility of ThunderRW. Experimental results show
that ThunderRW outperforms state-of-the-art approaches by an order of magnitude,
and the step interleaving technique significantly reduces the CPU pipeline stall
from 73.1% to 15.0%.

For the details, please refer to our VLDB'2021 paper
"ThunderRW: An In-Memory Graph Random Walk Engine"
by [Shixuan Sun](https://shixuansun.github.io/), [Yuhang Chen](https://alexcyh7.github.io/),
[Shengliang Lu](https://github.com/lushl9301), [Bingsheng He](https://www.comp.nus.edu.sg/~hebs/) and
[Yuchen Li](http://yuchenli.net/). If you have any further questions, please feel free to contact us.

Please cite our paper, if you use our source code.

* "Shixuan Sun, Yuhang Chen, Shengliang Lu, Bingsheng He and Yuchen Li.
ThunderRW: An In-Memory Graph Random Walk Engine. VLDB 2021."


## Compile

Under the root directory of the project, execute the following commands to compile the source code

```zsh
mkdir build
cd build
cmake ..
make
```

## Preprocessing

ThunderRW takes a graph stored in the CSR format as the input. It supports directed,
edge labeled and edge weighted graphs. Particularly, it requires that
the graph data is stored as follows: 1) the vertex ID is ranged from 0 to N-1
where N is the number of vertices in the graph; 2) the graph data such as
the graph structure, the edge label and the edge weight is stored in the same
folder; 3) b_degree.bin is an int32_t array recording the degree of each vertex where
the first element is unused, the second element records the number of vertices,
and the following elements are vertex degrees; 4) b_adj.bin is an int32_t array
recording the neighbors of each vertex (i.e., the neighbor array in CSR); 5)
b_edge_label is an int32_t array recording the label of each edge with the
one-to-one relationship to the b_adj.bin; and 6) b_edge_weight.bin is a double
array recording the weight of each edge with the one-to-one relationship to the
b_adj.bin. Therefore, you need to convert the graph file into the type that
can be consumed by ThunderRW.

We provide the script `prepare_data.sh` to convert the edge list file. Before
running our converting tool, you must change your graph file name
into `b_edge_list.bin`. In this repository, we use the amazon dataset as
the running example. The edge list file is under `sample_data/amazon` folder.
In the script, `dataset_path` sets the dataset root path, `array` stores
the dataset name list and `skip_char` configures the skip character. Execute
the following command to convert the data.

```zsh
./prepare_data.sh
```

You will see the following files in the data folder: `b_degree.bin`,
`b_adj.bin`, `b_edge_label.bin` and `b_edge_weight.bin`.
The tool converts the edge list to undirected graph and randomly assigned a label and a weight to each edge.
We use undirected graphs as the input to keep the workload of random walks queries nearly the same in our experiments.

## Execution

We implement four algorithms with ThunderRW, which include PPR, DeepWalk,
Node2Vec and MetaPath. The source files are under `random_walk/apps` folder.
The parameters have been configured in the source files. You can rewrite
the default setting (e.g., the number of walkers, the sampling
method and the execution mode) in the source files or through parameters with the
instructions in the files. Here, we demonstrate the input of graph data.

Execute PPR with the following command. `-f` sets the input graph folder and
`-n` sets the number of threads.

```zsh
./build/random_walk/ppr.out -f sample_dataset/amazon -n 10
```

Execute DeepWalk with the following command. `-ew` is to load the edge
weight array into main memory.

```zsh
./build/random_walk/deepwalk.out -f sample_dataset/amazon -n 10 -ew
```

Execute Node2Vec with the following command.

```zsh
./build/random_walk/node2vec.out -f sample_dataset/amazon -n 10 -ew
```

Execute MetaPath with the following command. `-el` is to load the
edge label array into main memory. `-s` is to set the meta path schema.

```zsh
./build/random_walk/metapath.out -f sample_dataset/amazon -n 10 -el -s 0,1,2,3,4
```

## Configuration

In `random_walk/types.h`, you can disable the step interleaving technique
by commenting out `ENABLE_STEPINTERLEAVING`. The ring size k and k´ can be configured
by `RING_SIZE` and `SEARCH_RING_SIZE`.

## Experiment Datasets

You can download the graphs used in our paper by following the
instructions in Section 6.1.