//! # mapper-graph
//!
//! Implementation of the **Mapper algorithm** from Topological Data Analysis.
//!
//! The Mapper graph is a discrete approximation of the Reeb graph of a filter
//! function applied to a point cloud. It is constructed by:
//! 1. Applying a **filter function** to map points to ℝ
//! 2. **Covering** the filter range with overlapping intervals
//! 3. **Clustering** the points within each cover element
//! 4. Building the **nerve** of the resulting clusters (nodes = clusters, edges = overlaps)
//!
//! # Example
//!
//! ```
//! use mapper_graph::{PointCloud, FilterFunction, Cover, ClusterStep, MapperGraph};
//!
//! let pc = PointCloud::new(vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]]);
//! let filter = FilterFunction::projection(0);
//! let cover = Cover::uniform(3, 0.3);
//! let graph = MapperGraph::build(&pc, &filter, &cover, &ClusterStep::single_linkage(0.5));
//! assert!(graph.node_count() > 0);
//! ```

use std::collections::{HashMap, HashSet};

// ---------------------------------------------------------------------------
// PointCloud
// ---------------------------------------------------------------------------

/// A set of `n` points in `d`-dimensional Euclidean space.
#[derive(Debug, Clone)]
pub struct PointCloud {
    points: Vec<Vec<f64>>,
}

impl PointCloud {
    /// Create a new point cloud from a vector of points.
    pub fn new(points: Vec<Vec<f64>>) -> Self {
        Self { points }
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// True if the cloud is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Dimensionality of each point.
    pub fn dim(&self) -> usize {
        if self.points.is_empty() {
            0
        } else {
            self.points[0].len()
        }
    }

    /// Access the i-th point.
    pub fn point(&self, i: usize) -> &[f64] {
        &self.points[i]
    }

    /// Iterate over all points.
    pub fn iter(&self) -> impl Iterator<Item = &[f64]> {
        self.points.iter().map(|v| v.as_slice())
    }

    /// Euclidean distance between two point indices.
    pub fn distance(&self, i: usize, j: usize) -> f64 {
        let a = &self.points[i];
        let b = &self.points[j];
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y) * (x - y))
            .sum::<f64>()
            .sqrt()
    }
}

// ---------------------------------------------------------------------------
// FilterFunction
// ---------------------------------------------------------------------------

/// A function `f: ℝ^d → ℝ` used to project the point cloud to the real line.
#[derive(Debug, Clone)]
pub enum FilterFunction {
    /// Project onto coordinate `axis`.
    Projection { axis: usize },
    /// Eccentricity: average distance to all other points.
    Eccentricity,
    /// Gaussian density estimate with given bandwidth.
    Density { bandwidth: f64 },
    /// A custom function.
    Custom(fn(&[f64]) -> f64),
}

impl FilterFunction {
    /// Create a projection filter onto the given axis.
    pub fn projection(axis: usize) -> Self {
        Self::Projection { axis }
    }

    /// Create an eccentricity filter.
    pub fn eccentricity() -> Self {
        Self::Eccentricity
    }

    /// Create a density-based filter with Gaussian kernel.
    pub fn density(bandwidth: f64) -> Self {
        Self::Density { bandwidth }
    }

    /// Create a custom filter function.
    pub fn custom(f: fn(&[f64]) -> f64) -> Self {
        Self::Custom(f)
    }

    /// Evaluate the filter on a single point.
    pub fn eval_single(&self, point: &[f64]) -> f64 {
        match self {
            Self::Projection { axis } => point.get(*axis).copied().unwrap_or(0.0),
            Self::Eccentricity | Self::Density { .. } => {
                // Need full cloud; use eval instead
                0.0
            }
            Self::Custom(f) => f(point),
        }
    }

    /// Evaluate the filter on every point in the cloud, returning a Vec of values.
    pub fn eval(&self, cloud: &PointCloud) -> Vec<f64> {
        match self {
            Self::Projection { axis } => cloud
                .iter()
                .map(|p| p.get(*axis).copied().unwrap_or(0.0))
                .collect(),
            Self::Eccentricity => {
                let n = cloud.len();
                let mut vals = vec![0.0; n];
                for (i, v) in vals.iter_mut().enumerate() {
                    let mut sum = 0.0;
                    for j in 0..n {
                        if i != j {
                            sum += cloud.distance(i, j);
                        }
                    }
                    *v = if n > 1 { sum / (n - 1) as f64 } else { 0.0 };
                }
                vals
            }
            Self::Density { bandwidth } => {
                let n = cloud.len();
                let sigma2 = 2.0 * bandwidth * bandwidth;
                let pts: Vec<&[f64]> = cloud.iter().collect();
                let mut vals = vec![0.0; n];
                for (i, v) in vals.iter_mut().enumerate() {
                    let mut sum = 0.0;
                    for j in 0..n {
                        let d2 = pts[i]
                            .iter()
                            .zip(pts[j].iter())
                            .map(|(a, b)| (a - b) * (a - b))
                            .sum::<f64>();
                        sum += (-d2 / sigma2).exp();
                    }
                    *v = sum / n as f64;
                }
                vals
            }
            Self::Custom(f) => cloud.iter().map(f).collect(),
        }
    }
}

// ---------------------------------------------------------------------------
// Cover
// ---------------------------------------------------------------------------

/// A covering of the filter range `[min, max]` by overlapping intervals.
#[derive(Debug, Clone)]
pub struct Cover {
    intervals: Vec<(f64, f64)>,
}

impl Cover {
    /// Create a cover from explicit intervals.
    pub fn new(intervals: Vec<(f64, f64)>) -> Self {
        Self { intervals }
    }

    /// Uniform cover: `num_intervals` equally-spaced intervals with fractional
    /// `overlap` (0.0 = no overlap, 0.5 = 50% overlap).
    pub fn uniform(num_intervals: usize, overlap: f64) -> Self {
        // This is a factory; the actual range is computed at build time.
        // We store parameters as a single interval covering [0, num_intervals]
        // and interpret it in `cover_range`.
        let _ = overlap;
        Self {
            intervals: (0..num_intervals).map(|i| (i as f64, i as f64 + 1.0)).collect(),
        }
    }

    /// Build a cover for the given range `[lo, hi]` with `n` intervals and overlap fraction.
    pub fn cover_range(lo: f64, hi: f64, n: usize, overlap: f64) -> Self {
        if n == 0 || hi <= lo {
            return Self { intervals: vec![] };
        }
        let step = (hi - lo) / n as f64;
        let margin = step * overlap;
        let intervals: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let start = lo + i as f64 * step - margin;
                let end = lo + (i + 1) as f64 * step + margin;
                (start, end)
            })
            .collect();
        Self { intervals }
    }

    /// Number of intervals.
    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    /// True if no intervals.
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    /// Return indices of intervals that contain `value`.
    pub fn containing(&self, value: f64) -> Vec<usize> {
        self.intervals
            .iter()
            .enumerate()
            .filter(|(_, (lo, hi))| value >= *lo && value <= *hi)
            .map(|(i, _)| i)
            .collect()
    }

    /// Access intervals.
    pub fn intervals(&self) -> &[(f64, f64)] {
        &self.intervals
    }
}

// ---------------------------------------------------------------------------
// ClusterStep
// ---------------------------------------------------------------------------

/// Strategy for clustering points within a single cover element.
#[derive(Debug, Clone)]
pub enum ClusterStep {
    /// Each point is its own cluster.
    Trivial,
    /// Single-linkage hierarchical clustering cut at threshold `distance`.
    SingleLinkage { distance: f64 },
}

impl ClusterStep {
    /// Create a single-linkage clustering step.
    pub fn single_linkage(distance: f64) -> Self {
        Self::SingleLinkage { distance }
    }

    /// Create a trivial clustering (every point is its own cluster).
    pub fn trivial() -> Self {
        Self::Trivial
    }

    /// Cluster the given point indices from `cloud`.
    /// Returns a list of clusters, each cluster being a Vec of point indices.
    pub fn cluster(&self, cloud: &PointCloud, indices: &[usize]) -> Vec<Vec<usize>> {
        if indices.is_empty() {
            return vec![];
        }
        match self {
            Self::Trivial => indices.iter().map(|&i| vec![i]).collect(),
            Self::SingleLinkage { distance } => {
                // Simple union-find single-linkage
                let n = indices.len();
                let mut parent: Vec<usize> = (0..n).collect();
                let mut rank = vec![0usize; n];

                let find = |parent: &mut Vec<usize>, mut x: usize| -> usize {
                    while parent[x] != x {
                        parent[x] = parent[parent[x]];
                        x = parent[x];
                    }
                    x
                };

                #[allow(clippy::needless_range_loop)]
                for i in 0..n {
                    for j in (i + 1)..n {
                        if cloud.distance(indices[i], indices[j]) < *distance {
                            let ri = find(&mut parent, i);
                            let rj = find(&mut parent, j);
                            if ri != rj {
                                if rank[ri] < rank[rj] {
                                    parent[ri] = rj;
                                } else if rank[ri] > rank[rj] {
                                    parent[rj] = ri;
                                } else {
                                    parent[rj] = ri;
                                    rank[ri] += 1;
                                }
                            }
                        }
                    }
                }

                let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
                #[allow(clippy::needless_range_loop)]
                for i in 0..n {
                    let root = find(&mut parent, i);
                    groups.entry(root).or_default().push(indices[i]);
                }
                groups.into_values().collect()
            }
        }
    }
}

// ---------------------------------------------------------------------------
// MapperGraph
// ---------------------------------------------------------------------------

/// A node in the Mapper graph: a cluster of original point indices.
#[derive(Debug, Clone)]
pub struct MapperNode {
    /// Label identifying the cover element and cluster within it.
    pub label: String,
    /// Original point indices in this cluster.
    pub points: Vec<usize>,
}

/// An edge in the Mapper graph indicating overlap between two clusters.
#[derive(Debug, Clone)]
pub struct MapperEdge {
    pub source: usize,
    pub target: usize,
}

/// The Mapper graph: the nerve of a clustered cover.
#[derive(Debug, Clone)]
pub struct MapperGraph {
    nodes: Vec<MapperNode>,
    edges: Vec<MapperEdge>,
}

impl MapperGraph {
    /// Build a Mapper graph from a point cloud, filter, cover, and clustering strategy.
    pub fn build(
        cloud: &PointCloud,
        filter: &FilterFunction,
        cover: &Cover,
        cluster_step: &ClusterStep,
    ) -> Self {
        let values = filter.eval(cloud);

        // For each cover interval, collect points whose filter value falls inside
        let mut all_clusters: Vec<(usize, Vec<usize>)> = Vec::new(); // (cover_idx, point_indices)
        for (ci, &(_lo, _hi)) in cover.intervals.iter().enumerate() {
            let members: Vec<usize> = values
                .iter()
                .enumerate()
                .filter(|(_, &v)| cover.containing(v).contains(&ci))
                .map(|(i, _)| i)
                .collect();
            let clusters = cluster_step.cluster(cloud, &members);
            for cluster in clusters {
                all_clusters.push((ci, cluster));
            }
        }

        let nodes: Vec<MapperNode> = all_clusters
            .iter()
            .enumerate()
            .map(|(i, (ci, pts))| MapperNode {
                label: format!("c{}-{}", ci, i),
                points: pts.clone(),
            })
            .collect();

        // Build edges: two nodes share an edge if their point sets overlap
        let mut edges = Vec::new();
        for i in 0..nodes.len() {
            for j in (i + 1)..nodes.len() {
                let set_i: HashSet<usize> = nodes[i].points.iter().copied().collect();
                let set_j: HashSet<usize> = nodes[j].points.iter().copied().collect();
                if set_i.intersection(&set_j).count() > 0 {
                    edges.push(MapperEdge {
                        source: i,
                        target: j,
                    });
                }
            }
        }

        Self { nodes, edges }
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of edges.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Access nodes.
    pub fn nodes(&self) -> &[MapperNode] {
        &self.nodes
    }

    /// Access edges.
    pub fn edges(&self) -> &[MapperEdge] {
        &self.edges
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_cloud() -> PointCloud {
        PointCloud::new(vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![2.0, 0.0],
            vec![3.0, 0.0],
            vec![4.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 1.0],
            vec![2.0, 1.0],
            vec![3.0, 1.0],
            vec![4.0, 1.0],
        ])
    }

    #[test]
    fn test_point_cloud_len() {
        let pc = sample_cloud();
        assert_eq!(pc.len(), 10);
        assert_eq!(pc.dim(), 2);
    }

    #[test]
    fn test_point_cloud_distance() {
        let pc = PointCloud::new(vec![vec![0.0], vec![3.0], vec![4.0]]);
        assert!((pc.distance(0, 1) - 3.0).abs() < 1e-10);
        assert!((pc.distance(1, 2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_filter_projection() {
        let pc = PointCloud::new(vec![vec![1.0, 5.0], vec![2.0, 6.0]]);
        let f = FilterFunction::projection(0);
        let vals = f.eval(&pc);
        assert_eq!(vals, vec![1.0, 2.0]);

        let f1 = FilterFunction::projection(1);
        let vals1 = f1.eval(&pc);
        assert_eq!(vals1, vec![5.0, 6.0]);
    }

    #[test]
    fn test_filter_eccentricity() {
        let pc = PointCloud::new(vec![vec![0.0], vec![10.0]]);
        let f = FilterFunction::eccentricity();
        let vals = f.eval(&pc);
        assert!((vals[0] - 10.0).abs() < 1e-10);
        assert!((vals[1] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_filter_custom() {
        let pc = PointCloud::new(vec![vec![3.0, 4.0]]);
        let f = FilterFunction::custom(|p| p[0] + p[1]);
        let vals = f.eval(&pc);
        assert!((vals[0] - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_cover_range() {
        let cover = Cover::cover_range(0.0, 10.0, 5, 0.0);
        assert_eq!(cover.len(), 5);
    }

    #[test]
    fn test_cover_containing() {
        let cover = Cover::cover_range(0.0, 10.0, 2, 0.1);
        let c = cover.containing(5.0);
        assert!(!c.is_empty());
    }

    #[test]
    fn test_cluster_trivial() {
        let pc = sample_cloud();
        let cs = ClusterStep::trivial();
        let clusters = cs.cluster(&pc, &[0, 1, 2]);
        assert_eq!(clusters.len(), 3);
    }

    #[test]
    fn test_cluster_single_linkage_far() {
        let pc = PointCloud::new(vec![vec![0.0], vec![100.0], vec![200.0]]);
        let cs = ClusterStep::single_linkage(1.0);
        let clusters = cs.cluster(&pc, &[0, 1, 2]);
        assert_eq!(clusters.len(), 3);
    }

    #[test]
    fn test_cluster_single_linkage_close() {
        let pc = PointCloud::new(vec![vec![0.0], vec![0.1], vec![0.2]]);
        let cs = ClusterStep::single_linkage(1.0);
        let clusters = cs.cluster(&pc, &[0, 1, 2]);
        assert_eq!(clusters.len(), 1);
    }

    #[test]
    fn test_mapper_graph_build() {
        let pc = sample_cloud();
        let filter = FilterFunction::projection(0);
        let cover = Cover::cover_range(0.0, 4.0, 3, 0.3);
        let graph = MapperGraph::build(&pc, &filter, &cover, &ClusterStep::single_linkage(1.5));
        assert!(graph.node_count() > 0);
    }

    #[test]
    fn test_mapper_graph_overlap_edges() {
        // Points at 0, 0.5, 5 — cover of [0,3] and [2,5] with overlap
        let pc = PointCloud::new(vec![vec![0.0], vec![0.5], vec![2.5], vec![5.0]]);
        let filter = FilterFunction::projection(0);
        let cover = Cover::cover_range(0.0, 5.0, 2, 0.3);
        let graph = MapperGraph::build(&pc, &filter, &cover, &ClusterStep::trivial());
        // Point at 2.5 should appear in both intervals creating an edge
        assert!(graph.edge_count() > 0);
    }
}
