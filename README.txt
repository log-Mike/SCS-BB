Size-bounded Community Search
This project aims to find a subgraph with the largest min-degree among all connected subgraphs that contain the query vertex q and have at least l and at most h vertices, where q, l, h are specified by the query. It outperforms the baseline algorithms by several orders of magnitude. Details can be found in our paper http://vldb.org/pvldb/vol14/p1441-yao.pdf

SCS is the executable, and is compiled on Ubuntu 18.04.5, with -O3 optimization.

Running Format
./SCS [1. graph name] [2. size LB: l] [3. size UB: h] [4. query vid: q] [5. Heuristic 1: 0/1] [6. Heuristic 2: 0/1] [7. Heuristic 3: 0/1] [8. UB1: 0/1] [9. UB2: 0/1] [10. UB3: 0/1] [11. UB3 optimization: 0/1] [12. reduction 1: 0/1] [13. reduction 2: 0/1] [14. reduction 3: 0/1] [15. dominating branch: 0/1] [16. dominating pair threshold: 2/4/8/16/...] [17. time limit: (in seconds)] [18. branching order: 1/2/3] 

Running Example
./SCS ./graph.txt 6 9 2 1 1 1 1 1 1 0 1 1 1 0 0 1800 1