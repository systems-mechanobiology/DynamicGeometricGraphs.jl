# DynamicGeometricGraphs.jl

A small package representing geometric graphs where vertices carry
N-dimensional coordinates and optional per-vertex metadata.

## Why This Package?

Built to address specific needs in geometric graph processing:

| Feature | This Package | MetaGraphs.jl | SimpleWeightedGraphs.jl |
|---------|--------------|---------------|-------------------------|
| N-dimensional coordinates | ✓ | ✗ | ✗ |
| Fast vertex/edge add/remove | ✓ | ~ | ✗ |
| Edit tracking | ✓ | ✗ | ✗ |
| Custom ambient metric | ✓ | ✗ | ✗ |
| On-the-fly edge weights | ✓ | ✗ | ✗ |
| Geometric transformations | ✓ | ✗ | ✗ |
| Sparse storage | ✓ | ✓ | ✗ |

**Best for:** Dynamic geometric graphs with frequent edits and coordinate-based distances.

**Not for:** Static graphs or arbitrary metadata schemas (use Graphs.jl or MetaGraphs.jl).

## Quick start

```julia
using StaticArrays, DynamicGeometricGraphs
using Graphs: add_vertex!, add_edge!

# 2D graph with Float64 coordinates and default Euclidean ambient metric
g = DynamicGeometricGraph{2, Float64}()

v1 = add_vertex!(g, SVector(0.0, 0.0))
v2 = add_vertex!(g, SVector(1.0, 0.0))
add_edge!(g, v1, v2)

# Edge weights are computed on-the-fly from coordinates
w = edge_weight(g, v1, v2)  # Returns 1.0

# get coordinates matrix and vertex ids
P, ids = get_coords(g)
```

## Ambient Metric

The **ambient metric** (distance function) determines how edge weights are computed from vertex coordinates:

- Set at construction via `distfun` parameter (default: Euclidean distance)
- Immutable once the graph is created
- Can be overridden per-call via `edge_weight(g, u, v; distancefun=custom_metric)`

```julia
# Custom metric at construction
using Graphs: add_vertex!, add_edge!

manhattan(a, b) = sum(abs.(a .- b))
g = DynamicGeometricGraph{2, Float64}(distfun=manhattan)

# Runtime override for specific queries
w = edge_weight(g, v1, v2; distancefun=euclid)
```

See the API page for exported functions and types.
