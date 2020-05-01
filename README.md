# Random Spanning Tree Sampler in C

## Build Instructions
``` bash
mkdir build && cd build
cmake ..
make
```
## Dependencies
- igraph: a C library for creating, manipulating and analysing graphs 
  - [GitHub](https://github.com/igraph/igraph)
  - [Installation Guide](https://igraph.org/c/)
- CMake and other build tools

## Currently Implemented Graphs
- Small graph generated by specifing edges
- Random connnected graph, for a given number of vertices, density target, maximum and minimum degree
- Full graph, Cycle graph, Lattice graph
- Many other graphs that are included in the igraph library [Graph Generators](https://igraph.org/c/doc/igraph-Generators.html)

## Limitations Found

- It is quite typical to have the sample generation time to be on the scale of 0.01ms. This will increase with the mean hitting time.
- The sheer amount of edges to calculate ratio with, when vertices increase, turns this algorithm not scalable.
- Sample reusable does help but has limited improvement over the scalability.
- IDEAS: could we estimation ratios of multiple edges together (with the same vertex, or far apart?), or establish some divide and conquer strategy so the hitting time for each sub-problem is dramatically decreased.

## Running Time of The Approximate Counter

### 2020-04-30 Update 2

- Bug fixes (invalid array index)
- Add: convergence test based on variance, instead of max/min

### 100 Vertex, Complete Graph

```
Thu Apr 30 18:55:40 2020
Created graph with 100 vertices and 4950 edges
1028731 samples taken, with per sample time taking 0.004 ms
Total time spent 4 seconds

FINAL result = 1.9114e+195 (e^4.4965e+02) with 6806777 effective samples, avg 1375 samples per edge, 
```

### 200 Vertex, Complete Graph

```
Thu Apr 30 18:58:07 2020
Created graph with 200 vertices and 19900 edges
3103882 samples taken, with per sample time taking 0.012 ms
Total time spent 36 seconds

FINAL result = inf (e^1.0499e+03) with 23226525 effective samples, avg 1167 samples per edge, 
```

### 500 Vertex, Complete Graph

```
Thu Apr 30 18:59:25 2020
Created graph with 500 vertices and 124750 edges
15563688 samples taken, with per sample time taking 0.044 ms
Total time spent 684 seconds

FINAL result = inf (e^3.0953e+03) with 115017822 effective samples, avg 921 samples per edge, 
```

true result = e^3094

### 1000 Vertex, Complete Graph

```
Thu Apr 30 19:32:09 2020
Created graph with 1000 vertices and 499500 edges
53717546 samples taken, with per sample time taking 0.118 ms
Total time spent 6335 seconds

FINAL result = inf (e^6.8905e+03) with 387723790 effective samples, avg 776 samples per edge. 
```

### 500 Vertex, 0.1 Density, Max Degree 5

```
Thu Apr 30 19:11:50 2020
Created graph with 500 vertices and 1218 edges
1281073 samples taken, with per sample time taking 0.039 ms
Total time spent 50 seconds

FINAL result = 2.5223e+305 (e^7.0321e+02) with 2760247 effective samples, avg 2266 samples per edge, 
```


### 1000 Vertex, 0.1 Density, Max Degree 5

```
Thu Apr 30 19:14:01 2020
Created graph with 1000 vertices and 2472 edges
2647359 samples taken, with per sample time taking 0.100 ms
Total time spent 265 seconds

FINAL result = inf (e^1.4328e+03) with 5701466 effective samples, avg 2306 samples per edge. 
```

observe that the mean hitting time become much worse

### 2020-04-30 Update

- Fix: avoid calculating ratio for edges removed by previous edge contraction. Previously causing slow down at later stage of the itereations
- Fix: Better estimation results now
- Improvement: constant factor speed boost

### 200 Vertex, Complete Graph

```
Thu Apr 30 00:01:02 2020
Created graph with 200 vertices and 19900 edges
3031779 samples taken, with per sample time taking 0.018 ms
Total time spent 55 seconds

FINAL result = inf (e^1.04e+03) with 17324838 effective samples
```

compared to the oldest implementation, 414 seconds.

Without randomising edges (edge sequenced according to vertex number)

```
Thu Apr 30 10:21:13 2020
Created graph with 200 vertices and 19900 edges
2950571 samples taken, with per sample time taking 0.011 ms
Total time spent 33 seconds

FINAL result = inf (e^1.05e+03) with 18427274 effective samples
```

### 500 Vertex, Complete Graph

```
Thu Apr 30 01:37:28 2020
Created graph with 500 vertices and 124750 edges
14831544 samples taken, with per sample time taking 0.063 ms
Total time spent 937 seconds

FINAL result = inf (e^3.06e+03) with 78708367 effective samples
```

Without randomising edges (edge sequenced according to vertex number)

```
Thu Apr 30 10:23:14 2020
Created graph with 500 vertices and 124750 edges
14817429 samples taken, with per sample time taking 0.037 ms
Total time spent 547 seconds

FINAL result = inf (e^3.07e+03) with 84095442 effective samples
```

### 200 Vertex, 0.1 Density, Max Degree 5

```
Thu Apr 30 00:52:47 2020
Created graph with 200 vertices and 473 edges
424509 samples taken, with per sample time taking 0.011 ms
Total time spent 4 seconds

FINAL result = 1.6136e+118 (e^2.72e+02) with 904954 effective samples
```

### 500 Vertex, 0.1 Density, Max Degree 5

```
Created graph with 500 vertices and 1218 edges
1070589 samples taken, with per sample time taking 0.052 ms
Total time spent 55 seconds

FINAL result = 3.0521e+305 (e^7.03e+02) with 2277627 effective samples
```

Without randomising edges (edge sequenced according to vertex number)

```
Thu Apr 30 10:34:23 2020
Created graph with 500 vertices and 1218 edges
1077820 samples taken, with per sample time taking 0.035 ms
Total time spent 37 seconds

FINAL result = 2.3998e+305 (e^7.03e+02) with 2308714 effective samples
```

### 1000 Vertex, 0.1 Density, Max Degree 5

Without randomising edges (edge sequenced according to vertex number)

```
Thu Apr 30 14:30:23 2020
Created graph with 1000 vertices and 2472 edges
2226236 samples taken, with per sample time taking 0.096 ms
Total time spent 213 seconds

FINAL result = inf (e^1.43e+03) with 4756279 effective samples

```

### 2020-04-25 Update

I have done the following major improvement
- Remove igraph dependency for heavy calculation part. Rewrite a small `GraphLite` class to handle the undirected graph better, more efficient
  - Stored as edge list and incident list of vertices
  - No hashing, linear storage; all operation used in the random walk could be done in O(1) vector entry access
  - When edge and vertex are removed, they are not physically deleted, so the index will not be messed up; instead they are properly labelled as deleted
  - Contracted edges make connecting vertices behaving as one, by reassigning edges. Previous attempt to keep the contracted vertices does not work well, as more time are wasted random walking on those existing edges
- igraph is still used for graph generation purpose
- Implemented Wilson's algorithm exactly, in the hope to help uniformity
- Fixed ratio estimating issue, when the edge we want to estimated is already contracted away (no longer exists). In this case, the ratio contributing is simply unity 1
- Enabled predefined edges shuffling to make progress of ratio estimation more even
  - TODO: change predefine edges when need, e.g. when non-estimated edges are being contracted, we can get their ratios for free!

### 200 Vertices, 0.1 Density, No Maximum Degree Limit
```
Sat Apr 25 20:39:29 2020
Created graph with 200 vertices and 19900 edges
3808316 samples taken, with per sample time taking 0.022 ms
Total time spent 84 seconds

FINAL result = inf with 47415646 effective samples
```
---

_Old Attempt Below..._

### 50 Vertices, Full Graph
```
Fri Apr 24 11:41:14 2020
Created graph with 50 vertices and 1225 edges
454403 samples taken, with per sample time taking 0.01 ms
Total time spent 6 seconds

FINAL result = 3.0410e+81 with 2541313 effective samples
```
### 100 Vertices, Full Graph
```
Fri Apr 24 11:26:45 2020
Created graph with 100 vertices and 4950 edges
1266882 samples taken, with per sample time taking 0.04 ms
Total time spent 53 seconds

FINAL result = 5.8982e+195 with 10333308 effective samples
```

### 200 Vertices, Full Graph
```
Fri Apr 24 11:42:09 2020
Created graph with 200 vertices and 19900 edges
3672212 samples taken, with per sample time taking 0.11 ms
Total time spent 414 seconds

FINAL result = inf with 41766822 effective samples
```
## Examples of Running Time of the Sampler

### 100 Vertices, 0.1 Density, No Maximum Degree Limit
```
Created graph with 100 vertices and 511 edges
...
100000 samples taken, with per sample time taking 36.17 us
Total time spent 3 seconds
```
This compares to the python implementation down from ~1ms per sample to 36us per sample. But still limited to linear performance improvement.

### 1000 Vertices
```
Created graph with 1000 vertices and 49714 edges
...
1000 samples taken, with per sample time taking 1.05 ms
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
```
Created graph with 10000 vertices and 49995000 edges
...
1000 samples taken, with per sample time taking 919.95 ms
Total time spent 919 seconds

```

MTT matrix method could do this in 90 seconds.

### Ring Graph, 1000 Vertices
```
Created graph with 1000 vertices and 1000 edges
...
1000 samples taken, with per sample time taking 16.95 ms
Total time spent 16 seconds
```

### 3D Lattice Graph, ～30000 Vertices
```
Created graph with 35937 vertices and 107811 edges
...
1000 samples taken, with per sample time taking 29.75 ms
Total time spent 29 seconds
```

Lattice seems to be friendlier to Markov Chains.