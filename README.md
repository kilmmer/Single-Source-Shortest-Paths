# Single-Source Shortest Paths (SSSP) Algorithm

This repository contains an implementation of the Single-Source Shortest Paths (SSSP) algorithm as described in the paper "Breaking the Sorting Barrier for Directed Single-Source Shortest Paths" by Ran Duan, Jiayi Mao, Xiao Mao, Xinkai Shu, and Longhui Yin (arXiv:2504.17033v2, 30 Jul 2025).

## Overview
The SSSP algorithm computes the shortest paths from a single source vertex to all other vertices in a directed graph with non-negative real weights. This implementation is written in TypeScript and is designed to run in Node.js environments.

### Key Features
- **Graph Representation**: The graph is represented as a `Map<number, {to: number, w: number}[]>`, where each key is a vertex and the value is an array of edges with their respective weights.
- **Recursive BMSSP Function**: Implements the BMSSP algorithm described in the paper, which uses boundary reduction and partial sorting.
- **Priority Queue**: A custom priority queue with `decreaseKey` functionality for efficient base case handling.
- **Partial Sorting Data Structure**: A block-based data structure for partial sorting, used in the recursive BMSSP function.
- **Logging**: Includes detailed logging functionality for debugging, with options to limit logs and disable logging.

## Usage
### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/kilmmer/Single-Source-Shortest-Paths
   ```
2. No dependencies are required for this project.

### Running the Algorithm
1. Define a graph:
   ```typescript
   const graph: Graph = {
     n: 4,
     adj: new Map([
       [0, [{to: 1, w: 1}, {to: 2, w: 4}]],
       [1, [{to: 2, w: 2}, {to: 3, w: 5}]],
       [2, [{to: 3, w: 1}]],
     ])
   };
   ```
2. Run the SSSP algorithm:
   ```typescript
   const distances = sssp(graph, 0);
   console.log(distances); // Expected output: [0, 1, 3, 4]
   ```

### Disabling Logs
To disable logs, run the script with the `--no-log` flag:
```bash
node sssp.js --no-log
```

## Limitations
- **Non-Negative Weights**: The algorithm assumes all edge weights are non-negative.
- **No Ties in Distances**: Assumes no ties in numerical distances, as per the paper.
- **Depth Limitations**: Recursive implementation may encounter stack depth limitations for very large graphs.

## References
- [Breaking the Sorting Barrier for Directed Single-Source Shortest Paths](https://arxiv.org/abs/2504.17033v2)

## Next Steps
I will continue studying, understanding, and refining the algorithm.

## License
This project is licensed under the MIT License.
