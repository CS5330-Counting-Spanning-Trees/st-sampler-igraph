# Random Spanning Tree Sampler in C
## Dependencies
- igraph: a C library for creating, manipulating and analysing graphs 
  - [GitHub](https://github.com/igraph/igraph)
  - [Installation Guide](https://igraph.org/c/)
- CMake and other build tools

## Currently Implemented Graphs
- Small graph generated by specifing edges
- Random connnected graph, for a given number of vertices, density target, maximum and minimum degree
- Full graph
- Many other graphs that are included in the igraph library [Graph Generators](https://igraph.org/c/doc/igraph-Generators.html)


## Examples of Running Time
### 1000 Vertices
```
Created graph with 1000 vertices and 49748 edges
...
1000 samples taken, with per sample time taking 1029.00 us
Total time spent 1 seconds
```
### 10000 Vertices
```
Created graph with 10000 vertices and 24974 edges
...
1000 samples taken, with per sample time taking 18.21 ms
Total time spent 18 seconds
```

### Complete Graph, 10000 Vertices