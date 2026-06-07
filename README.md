# mapper-graph

> **The Mapper algorithm from topological data analysis — point clouds → simplicial complexes that reveal shape**

[![crates.io](https://img.shields.io/crates/v/mapper-graph.svg)](https://crates.io/crates/mapper-graph)
[![docs.rs](https://docs.rs/mapper-graph/badge.svg)](https://docs.rs/mapper-graph)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## What is the Mapper Algorithm?

Imagine you have a cloud of points — sensor readings, gene expression data, customer purchase vectors — and you want to understand its **shape**. Are there clusters? Loops? Branches? The Mapper algorithm is a technique from **Topological Data Analysis (TDA)** that constructs a combinatorial graph (a simplicial complex) from any point cloud, revealing the underlying topological structure without assumptions about the data distribution.

Introduced by Singh, Mémoli, and Carlsson in 2007, the Mapper algorithm works by:
1. Projecting your high-dimensional data onto a single axis (the "filter function")
2. Covering that axis with overlapping intervals
3. Clustering the points within each interval
4. Connecting clusters that share points

The result is a graph whose topology mirrors the topology of your data — loops in the data appear as loops in the Mapper graph, clusters appear as dense regions, and branches appear as bifurcations.

## Why Does This Matter?

The Mapper algorithm has been applied to discover new breast cancer subtypes (Nicolau et al., 2011), analyze the shape of natural image patches (Carlsson et al., 2008), and study viral evolution trajectories. It's powerful because:

- **Model-free**: No assumptions about Gaussian distributions or linearity
- **Multi-scale**: The overlapping cover naturally captures features at different resolutions
- **Visual**: The resulting graph is immediately interpretable
- **Robust**: Small perturbations in the input produce small changes in the output (stability theorem)

Real-world applications include:
- **Biology**: Identifying disease subtypes from patient genomic data
- **Machine learning**: Feature engineering via topological descriptors
- **Data visualization**: Interactive exploration of high-dimensional datasets
- **Anomaly detection**: Unexpected topological features indicate outliers or novel patterns

## Architecture

```
┌──────────────────────────────────────────────────────────────┐
│                     Mapper Algorithm Pipeline                 │
│                                                              │
│  Point Cloud          Filter Function        Cover           │
│  ┌───┐    ┌──┐       ┌──────────┐     ┌───────────────┐    │
│  │ · │    │··│       │ f: ℝᵈ→ℝ  │     │ [a₁,b₁]      │    │
│  │· ·│ ──▶│··│ ────▶ │          │ ──▶ │ [a₂,b₂]      │    │
│  │ · │    │··│       │ e.g.     │     │ [a₃,b₃]      │    │
│  └───┘    └──┘       │ project  │     │  (overlap)    │    │
│                        │ density  │     └──────┬────────┘    │
│  Input data            │ custom   │            │             │
│                        └──────────┘            ▼             │
│                                            Cluster          │
│  ┌─────────────────────┐            ┌───────────────┐       │
│  │   Mapper Graph      │            │ SingleLinkage │       │
│  │   ┌───┐   ┌───┐    │ ◀────────  │ Trivial       │       │
│  │   │ C₁│───│ C₂│    │   Nerve    └───────────────┘       │
│  │   └─┬─┘   └───┘    │                                    │
│  │     │               │   Nodes = Clusters                 │
│  │   ┌─▼─┐            │   Edges = Overlapping clusters      │
│  │   │ C₃│            │                                    │
│  │   └───┘            │                                    │
│  └─────────────────────┘                                    │
└──────────────────────────────────────────────────────────────┘
```

## Quick Start

```rust
use mapper_graph::{PointCloud, FilterFunction, Cover, ClusterStep, MapperGraph};

// Create a point cloud (e.g., points along a circle)
let pc = PointCloud::new(vec![
    vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0],
    vec![1.0, 1.0], vec![0.5, 0.5],
]);

// Choose a filter function — here, project onto the x-axis
let filter = FilterFunction::projection(0);

// Cover the filter range with 3 overlapping intervals (30% overlap)
let cover = Cover::cover_range(0.0, 1.0, 3, 0.3);

// Cluster within each interval using single-linkage at distance 0.5
let cluster = ClusterStep::single_linkage(0.5);

// Build the Mapper graph
let graph = MapperGraph::build(&pc, &filter, &cover, &cluster);

println!("Nodes: {} (clusters of points)", graph.node_count());
println!("Edges: {} (overlapping clusters)", graph.edge_count());

// Inspect individual nodes
for node in graph.nodes() {
    println!("  {} contains {} points", node.label, node.points.len());
}
```

### Using Different Filter Functions

```rust
// Eccentricity: average distance to all other points
let filter = FilterFunction::eccentricity();

// Gaussian density estimate
let filter = FilterFunction::density(1.0);

// Custom function: distance from origin
let filter = FilterFunction::custom(|p| (p[0].powi(2) + p[1].powi(2)).sqrt());
```

## API Reference

### Point Cloud

| Method | Returns | Description |
|--------|---------|-------------|
| `PointCloud::new(points)` | `PointCloud` | Create a point cloud from `Vec<Vec<f64>>` |
| `pc.len()` | `usize` | Number of points |
| `pc.dim()` | `usize` | Dimensionality of each point |
| `pc.point(i)` | `&[f64]` | Access the i-th point |
| `pc.distance(i, j)` | `f64` | Euclidean distance between points i and j |

### FilterFunction

| Variant | Description |
|---------|-------------|
| `Projection { axis }` | Project onto coordinate axis |
| `Eccentricity` | Average distance to all other points |
| `Density { bandwidth }` | Gaussian kernel density estimate |
| `Custom(fn(&[f64]) -> f64)` | User-defined function |

### Cover

| Method | Description |
|--------|-------------|
| `Cover::cover_range(lo, hi, n, overlap)` | `n` intervals over `[lo, hi]` with fractional overlap |
| `cover.containing(value)` | Indices of intervals containing the value |

### ClusterStep

| Variant | Description |
|---------|-------------|
| `Trivial` | Each point is its own cluster |
| `SingleLinkage { distance }` | Hierarchical single-linkage cut at threshold |

### MapperGraph

| Method | Returns | Description |
|--------|---------|-------------|
| `MapperGraph::build(cloud, filter, cover, cluster)` | `MapperGraph` | Run the full Mapper pipeline |
| `graph.node_count()` | `usize` | Number of clusters (nodes) |
| `graph.edge_count()` | `usize` | Number of overlaps (edges) |
| `graph.nodes()` | `&[MapperNode]` | Access all nodes (each has `.label` and `.points`) |
| `graph.edges()` | `&[MapperEdge]` | Access all edges (source → target) |

## Mathematical Background

The Mapper algorithm approximates the **Reeb graph** of a function f: X → ℝ. Given:
- A topological space X (your point cloud)
- A filter function f: X → ℝ
- A cover 𝒰 = {U₁, ..., Uₙ} of f(X) with open sets

The Mapper graph M(f, 𝒰) is the nerve of the pullback cover:
```
M(f, 𝒰) = Nerve({f⁻¹(Uᵢ) : Uᵢ ∈ 𝒰})
```

After clustering each pullback f⁻¹(Uᵢ) into clusters Cᵢ₁, Cᵢ₂, ..., the nerve construction creates:
- A **node** for each cluster Cᵢⱼ
- An **edge** between Cᵢⱼ and Cₖₗ if they share at least one original data point

The **stability theorem** (Carrière & Oudot, 2018) guarantees that small changes in the input produce bounded changes in the Mapper output, measured by the bottleneck distance on the corresponding Reeb graphs.

### Parameters That Matter

| Parameter | Effect | Typical Range |
|-----------|--------|---------------|
| Number of intervals | Resolution of the graph | 5–50 |
| Overlap fraction | Connectivity between regions | 0.1–0.5 |
| Clustering threshold | Granularity within regions | Data-dependent |
| Filter function | Which aspect of shape to capture | Domain-dependent |

## Installation

```bash
cargo add mapper-graph
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
mapper-graph = "0.1.0"
```

## Related Crates

- [`cech-complex`](https://github.com/SuperInstance/cech-complex) — Čech complex via ball intersection nerve
- [`witness-complex`](https://github.com/SuperInstance/witness-complex) — Landmark-based topological approximation
- [`persistence-landscape`](https://github.com/SuperInstance/persistence-landscape) — Persistence landscapes for statistical TDA
- [`betti-curve`](https://github.com/SuperInstance/betti-curve) — Betti curves and Euler characteristic curves

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project — persistent cognitive substrate for multi-agent systems.*
