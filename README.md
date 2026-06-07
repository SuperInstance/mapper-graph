# mapper-graph

> **The Mapper algorithm from topological data analysis**

[![crates.io](https://img.shields.io/crates/v/mapper-graph.svg)](https://crates.io/crates/mapper-graph)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Implements the Mapper algorithm for constructing simplicial complexes from point clouds. Given a filter function, a covering of the filter range, and a clustering step, Mapper builds a graph that captures the topological structure of the data.

## Algorithm

1. **Filter**: Apply a function f: X → ℝ to the point cloud
2. **Cover**: Cover the range f(X) with overlapping intervals
3. **Cluster**: Within each interval, cluster the preimage points
4. **Nerve**: Build a graph where nodes = clusters, edges = overlapping clusters

The resulting Mapper graph is a topological summary that preserves features like clusters, loops, and branches.

## Installation

```toml
[dependencies]
mapper-graph = "0.1.0"
```

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project — persistent cognitive substrate for multi-agent systems.*
