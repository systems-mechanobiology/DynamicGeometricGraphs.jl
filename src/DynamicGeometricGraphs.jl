# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Ben Cardoen <b.cardoen@bham.ac.uk>

module DynamicGeometricGraphs

"""
DynamicGeometricGraphs

A Julia package representing geometric graphs where vertices
have N-dimensional coordinates (stored as `SVector`) and edges are created
and weighted according to an ambient metric. DGG's are designed to have streams of edits, 
to support this it keeps track of edit count (resets on copy). 
For performance it leverages StaticArrays so operations on coordinates are in principle
allocation free and can be stack only.

Primary features
- a mutable graph type `DynamicGeometricGraph` storing coordinates and
    per-vertex metadata
- adjacency stored in nested dictionaries for fast lookup and mutation
- pluggable ambient metric (default: `euclid`)
- edge weights computed on-the-fly (no storage overhead)
- N-dimensional coordinates, with metadata.

This module implements the interface expected by `Graphs.jl`
so `DynamicGeometricGraph` can be used with Graphs utilities.

Example
```
using StaticArrays
g = DynamicGeometricGraph{2, Float64}()
v = add_vertex!(g, SVector(0.0, 0.0))
```

Note: In the ecosystem there are SpatialGraphs, and MetaGraphs, but these do not quite match our use case
of a fast, fixed sized, n-dimensional graph type designed for frequent edits.

"""

using Graphs: AbstractGraph, nv
using Graphs
using StaticArrays
using Bijections
using LinearAlgebra
using SparseArrays

export DynamicGeometricGraph,
       edge_weight,
       euclid,
       find_nearest,
       freeze,
       generate_hexagonal_points,
       generate_hub_spoke_graph,
       get_coords,
       get_metadata,
       get_neighbour_coords,
       get_unused_idx!,
       get_vertex_coords,
       get_vertex_idx,
       graph_edits,
       geomfilter,
       has_vertex,
       incident_edges,
       merge_graphs,
       precision,
       refgraph,
       rem_edge!,
       rotate_graph,
       scale_graph,
       shape_distance,
       similarity,
       set_metadata!,
       sort_clockwise_indices,
       translate_graph,
       update_coord!,
       update_graph,
       weighted_degree

"""
DynamicGeometricGraph{N, T, F}

Type representing a dynamic geometric graph with N-dimensional vertex
coordinates of numeric type `T` and an ambient metric of type `F`.

Fields
- `vertices::Dict{Int, SVector{N,T}}`: map vertex id -> coordinates
- `metadata::Dict{Int, SVector{N,T}}`: optional per-vertex metadata
- `lookupvertices::Dict{SVector{N,T}, Int}`: reverse map coords -> id
- `edges::Dict{Int, Set{Int}}`: adjacency lists (weights computed on-the-fly)
- `dim::Int`: dimensionality N
- `counter::MVector{2,Int}`: internal counters (next id, edit count)
- `precision::T`: numeric precision used for metadata defaults
- `distance::F`: ambient metric function used to compute spatial distance/weight

This type implements `Graphs.AbstractGraph` interface so it can
be used with functions from the Graphs.jl ecosystem.

Edge weights are computed on-the-fly from vertex coordinates using the
ambient metric - no weight storage overhead.
"""
struct DynamicGeometricGraph{N, T<:Real, F<:Function} <: Graphs.AbstractGraph{Int}
    vertices :: Dict{Int, SVector{N, T}}       # store each vertex's coordinates
    metadata :: Dict{Int, SVector{N, T}}
    lookupvertices :: Dict{SVector{N, T}, Int} # Reverse map
    edges    :: Dict{Int, Set{Int}}            # adjacency info: neighbors only
    # forest   :: Vector{NearestNeighbors.BallTree}
    dim      :: Int
    counter  :: MVector{2,Int}
    precision :: T
    distance :: F

    # Primary constructor
    function DynamicGeometricGraph{N, T, F}(vertices::Dict{Int, SVector{N, T}},
        metadata::Dict{Int, SVector{N, T}},
        lookupvertices::Dict{SVector{N, T}, Int},
        edges::Dict{Int, Set{Int}},
        # forest::Vector{BallTree},
        dim::Int,
        counter::MVector{2,Int},
        precision ::T,
        distfun::F) where {N, T<:Real, F<:Function}
    return new(vertices, metadata, lookupvertices, edges, dim, counter, precision, distfun)
    end
end

function DynamicGeometricGraph{N, T}(;distfun = euclid, precision=1.0) where {N, T<:Real}
    DynamicGeometricGraph{N,T,typeof(distfun)}(
        Dict{Int,SVector{N,T}}(),
        Dict{Int,SVector{N,T}}(),
        Dict{SVector{N,T}, Int}(),
        Dict{Int, Set{Int}}(),
        N,
        MVector{2,Int}(0, 0),
        precision,
        distfun
    )
end


function incident_edges(g::DynamicGeometricGraph{N, T, F}, u::Int) where {N, T, F}
    if !haskey(g.edges, u)
        return []
    else
        return ((u, v, edge_weight(g, u, v)) for v in g.edges[u])
    end
end

Graphs.outneighbors(g::DynamicGeometricGraph, u::Int) = _nbr(g, u)


Graphs.inneighbors(g::DynamicGeometricGraph, u::Int) = _nbr(g, u)

function _nbr(g::DynamicGeometricGraph, u::Int)
    (v for (_, v, _) in incident_edges(g, u))
end

# Implement required Graphs.jl interface methods
Graphs.nv(g::DynamicGeometricGraph) = length(g.vertices)
Graphs.ne(g::DynamicGeometricGraph) = isempty(g.edges) ? 0 : sum(length(neighbors) for neighbors in values(g.edges)) ÷ 2
Graphs.vertices(g::DynamicGeometricGraph) = keys(g.vertices)
Graphs.edges(g::DynamicGeometricGraph) = (Graphs.Edge(u, v) for u in keys(g.edges) for v in g.edges[u] if u < v)
Graphs.is_directed(::Type{<:DynamicGeometricGraph}) = false
Graphs.is_directed(::DynamicGeometricGraph) = false

"""
edge_weight(g, u, v; adjacency=:spatial, distancefun=nothing)

Return the weight of the undirected edge between vertices `u` and `v`.

Edge weights are computed on-the-fly from vertex coordinates or metadata.

Parameters
- `adjacency`: `:spatial` (default) computes weight from vertex coordinates
  using `g.distance` (the graph's ambient metric). `:shape` computes weight 
  from `metadata` using `shape_distance` or an optional `distancefun`.
  **Note:** `:shape` mode is experimental and its API may change.
- `distancefun`: optional function to override the graph's ambient metric for this call.

Examples
```julia
# Use the graph's ambient metric
w = edge_weight(g, u, v)

# Override with custom metric for this call
manhattan(a, b) = sum(abs.(a .- b))
w = edge_weight(g, u, v; distancefun=manhattan)

# Use shape-based distance from metadata (EXPERIMENTAL)
w = edge_weight(g, u, v; adjacency=:shape)
```

Note: The ambient metric set at graph construction is immutable. Use `distancefun` 
to override the metric for specific edge weight queries.

!!! warning "Experimental Feature"
    The `adjacency=:shape` mode is experimental and under development. 
    The API may change in future versions.
"""
function edge_weight(g::DynamicGeometricGraph{N, T, F}, u::Int, v::Int; adjacency=:spatial, distancefun=nothing) where {N, T, F}
    if !Graphs.has_edge(g, u, v)
        @warn "Not a real edge"
    end
    if adjacency == :spatial
        return _edge_weight_spatial(g, u, v; distancefun=distancefun)
    end
    if adjacency == :shape 
        return _edge_weight_shape(g, u, v; distancefun=distancefun)
    end
    @assert false
end

function _edge_weight_shape(g::DynamicGeometricGraph{N, T, F}, u::Int, v::Int; distancefun=nothing) where {N, T, F}
    @warn "Shape-based adjacency is experimental. API may change in future versions." maxlog=1
    ud, vd = g.metadata[u], g.metadata[v]
    if isnothing(distancefun)
        return shape_distance(ud, vd; σ=g.precision)
    else
        return distancefun(ud, vd)
    end
end

function _edge_weight_spatial(g::DynamicGeometricGraph{N, T, F}, u::Int, v::Int; distancefun=nothing) where {N, T, F}
    if isnothing(distancefun)
        return g.distance(g.vertices[u], g.vertices[v])
    else
        return distancefun(g.vertices[u], g.vertices[v])
    end
end


Graphs.has_vertex(g::DynamicGeometricGraph, u) = DynamicGeometricGraphs._has_vertex(g, u)

function _has_vertex(g::DynamicGeometricGraph{N, T, F}, u::Int) where {N, T, F}
    haskey(g.vertices, u)
end

Graphs.has_edge(g::DynamicGeometricGraph, u, v) = DynamicGeometricGraphs._has_edge(g, u, v)

function _has_edge(g::DynamicGeometricGraph{N, T, F}, u::Int, v::Int) where {N, T, F}
    u ∈ keys(g.edges) && v ∈ g.edges[u]
end

function has_vertex(g::DynamicGeometricGraph{N, T, F}, v::SVector{N, T}) where {N, T, F}
    haskey(g.lookupvertices, v)
end

function get_vertex_coords(g::DynamicGeometricGraph{N, T, F}, u::Int) where {N, T, F}
    return get(g.vertices, u, nothing)
end

function get_vertex_idx(g::DynamicGeometricGraph{N, T, F}, v::SVector{N, T}) where {N, T, F}
    return get(g.lookupvertices, v, nothing)
end

function get_metadata(g::DynamicGeometricGraph, vertex::Int)
    get(g.metadata, vertex, nothing)
end

function set_metadata!(g::DynamicGeometricGraph, vertex::Int, meta)
    if haskey(g.metadata, vertex)
        g.metadata[vertex] = meta
        update_graph(g)
    else
        @warn "No such key $vertex, ignoring"
    end
end


"""
similarity(a, b; σ = exp(-1/2))

Simple similarity kernel used by shape-distance computations. Values close
to 1 indicate high similarity; values near 0 indicate dissimilarity.
"""
function similarity(a, b;σ =exp(-1/2))
    exp(- abs(a-b)/σ)
end

"""
shape_distance(a, b; σ=exp(-1/2))

Compute a shape-based distance in [0,1] using a similarity kernel.
`shape_distance` wraps `similarity` and returns `1 - similarity`. Both
scalar and vector metadata variants are supported. For vector metadata the
first element is used (with a warning).

!!! warning "Experimental"
    This function is part of the experimental shape-based adjacency feature.
    API may change in future versions.
"""
function shape_distance(a::T, b::T; σ=exp(-1/2)) where T<:Real
    1-similarity(a, b; σ=σ)
end


function shape_distance(a::AbstractVector{T}, b::AbstractVector{T}; σ=exp(-1/2)) where T<:Real
    @warn "Vector metadata, using first entry"
    1-similarity(a[1], b[1]; σ=σ)
end

"""
freeze(g)

Return a deep copy of `g` that can be mutated independently of the
original. Useful for creating a snapshot of the graph state.
"""
function freeze(g0::DynamicGeometricGraph{N, T, F}) where {N, T, F}
    sg = DynamicGeometricGraph{N, T, typeof(g0.distance)}(
        deepcopy(g0.vertices),
        deepcopy(g0.metadata),
        deepcopy(g0.lookupvertices),
        deepcopy(g0.edges),
        N,
        MVector{2,Int}(0, 0),
        g0.precision,
        g0.distance
    )
    sg.counter[1]=g0.counter[1]
    return sg
end

Graphs.rem_vertex!(g::DynamicGeometricGraph, u) = DynamicGeometricGraphs._rem_vertex!(g, u)

function _rem_vertex!(g::DynamicGeometricGraph{N, T, F}, u::Int) where {N, T, F}
    if Graphs.has_vertex(g, u)
        c = g.vertices[u]
        ies = incident_edges(g, u) |> collect
        s = pop!(g.lookupvertices, c, nothing)
        r = pop!(g.metadata, u, nothing)
        t = pop!(g.vertices, u, nothing)
        if isnothing(t) || isnothing(s) || isnothing(r)
            @error "Invalid vertex"
            return 
        else
            # ies = incident_edges(g, u)
            @debug "Dropping edges for $(u)"
            c=0
            for (u, v, _) in ies
                pop!(g.edges[u], v, nothing)
                pop!(g.edges[v], u, nothing)
                @debug "Dropping edge entry for $(v)"
                c+=1
            end 
            pop!(g.edges, u, nothing)
            g.counter[2] += c+1
        end
    else
        @warn "Deleting vertex that does not exist $(u)"
    end
end

"""
    Graphs.rem_edge!(g::DynamicGeometricGraph, u, v)

Remove the undirected edge between vertices `u` and `v` from graph `g`.

If the edge does not exist, a debug message is logged. If either vertex
does not exist, a warning is logged and no action is taken.

# Arguments
- `g::DynamicGeometricGraph`: The graph to modify
- `u::Int`: Source vertex ID
- `v::Int`: Destination vertex ID

# Returns
Nothing. Modifies `g` in-place.

# Example
```julia
g = DynamicGeometricGraph{2, Float64}()
v1 = add_vertex!(g, SVector(0.0, 0.0))
v2 = add_vertex!(g, SVector(1.0, 0.0))
add_edge!(g, v1, v2)
rem_edge!(g, v1, v2)  # Remove the edge
```
"""
Graphs.rem_edge!(g::DynamicGeometricGraph, u, v) = DynamicGeometricGraphs._rem_edge!(g, u, v)

function _rem_edge!(g::DynamicGeometricGraph{N, T, F}, u::Int, v::Int) where {N, T, F}
    if !Graphs.has_vertex(g, u) || !Graphs.has_vertex(g, v)
        @debug "$u or $v do not exist in g"
        return
    end
    if !haskey(g.edges, u) || !(v in g.edges[u])
        @debug "Edge $u -> $v did not exist"
        return
    end
    delete!(g.edges[u], v)
    haskey(g.edges, v) && delete!(g.edges[v], u)
    g.counter[2] += 1
end


"""
    Graphs.adjacency_matrix(g::DynamicGeometricGraph; adjacency=:spatial, distancefun=nothing)

Return the weighted adjacency matrix for a DynamicGeometricGraph.

This function builds a symmetric sparse matrix where entry (i,j) contains the edge weight
between vertices i and j. The matrix is constructed using the upper triangle only and then
made symmetric.

# Arguments
- `g::DynamicGeometricGraph`: The graph to extract the adjacency matrix from
- `adjacency::Symbol`: Type of adjacency to use (`:spatial` or `:shape`), default `:spatial`
- `distancefun::Function`: Optional custom distance function

# Returns
- A `LinearAlgebra.Symmetric` sparse matrix of edge weights

# Example
```julia
using DynamicGeometricGraphs
using Graphs

g = refgraph(100.0, 30)
A = adjacency_matrix(g)  # Weighted adjacency matrix
```
"""
function Graphs.adjacency_matrix(g::DynamicGeometricGraph{N, T, F}; adjacency=:spatial, distancefun=nothing) where {N, T, F}
    I = Int[]
    J = Int[]
    V = T[]
    n = nv(g)

    if n == 0
        @error "Empty graph"
        return sparse(I, J, V, 0, 0)
    end

    # Create vertex bijection: map vertex IDs to matrix indices [1, 2, ..., n]
    vs = sort(collect(Graphs.vertices(g)))
    vertex_to_idx = Dict(vid => idx for (idx, vid) in enumerate(vs))

    # Iterate through edges and collect (i, j, weight) triplets
    for (u_vid, neighbors) in g.edges
        if haskey(vertex_to_idx, u_vid)
            u_idx = vertex_to_idx[u_vid]
            for v_vid in neighbors
                if haskey(vertex_to_idx, v_vid)
                    v_idx = vertex_to_idx[v_vid]
                    # Only store upper triangle to avoid duplicates
                    if u_idx <= v_idx
                        push!(I, u_idx)
                        push!(J, v_idx)
                        w = edge_weight(g, u_vid, v_vid; adjacency=adjacency, distancefun=distancefun)
                        push!(V, w)
                    end
                end
            end
        end
    end

    # Build sparse matrix and make it symmetric
    A = sparse(I, J, V, n, n)
    A = Symmetric(A, :U)

    return A
end




"""
    Graphs.add_vertex!(g, coords)

Add a vertex using an N-dimensional coordinate vector coords.
Adding an existing vertex simply retrieves the index of the current vertex.
"""
function Graphs.add_vertex!(g::DynamicGeometricGraph{N,T}, coords::SVector{N, T};meta::SVector{N, T}=SVector{N, T}(Inf for _ in 1:N)) where {N, T<:Real}
    # @warn "Adding $coords to g"
    idx = get(g.lookupvertices, coords, -1)
    if any(isinf.(meta))
        # @info "Setting to precision"
        meta = SVector{N, T}(g.precision for _ in 1:N)
    end
    if idx == -1
        # new point, so add it with new index
        idx = get_unused_idx!(g)
        @debug "Adding $coords to $idx"
        g.lookupvertices[coords]=idx
        g.vertices[idx]=coords
        g.metadata[idx]=meta
        update_graph(g)
        # g.counter[2] += 1
    else
        @debug "Point exists, returning index $idx"
    end
    return idx
end

"""
    get_unused_idx!(g)
    
    For internal use, get a new index so during removal/adding of vertices we are always sure we do not use potentially invalidated vertex ids.
"""
function get_unused_idx!(g::DynamicGeometricGraph)
    g.counter[1]+=1
    return g.counter[1]
end

"""
    Graphs.add_vertex!
"""
function Graphs.add_vertex!(g::DynamicGeometricGraph{N,T}, coords::SVector{N, T}, index::Int) where {N, T<:Real}
    @error "Unsupported method"
    error(-1)
    # If I support this, I need to update all edges involving index
    if haskey(g.vertices, index)
        #this will invalidate edges
        error(-1)
    else
        if index < g.counter[1]
            # reclaiming index, but can be true that coords exists for another vertex
            if !haskey(g.lookupvertices, coords)
                g.vertices[index] = coords
                g.lookupvertices[coords] = index
            else 
                error(-1)
            end
            # Reclaiming an old index
        else
            error(-1)
        end
    end
end

"""
    Graphs.add_edge!(g, u, v)

Add an undirected edge between existing vertices u and v.
Edge weights are computed on-the-fly via the graph's ambient metric.
"""
function Graphs.add_edge!(g::DynamicGeometricGraph{N,T}, u::Int, v::Int) where {N, T<:Real}
    @assert haskey(g.vertices, u) && haskey(g.vertices, v)
    eset = get!(Set{Int}, g.edges, u)
    push!(eset, v)
    eset = get!(Set{Int}, g.edges, v)
    push!(eset, u)
    update_graph(g)
end


"""
precision(g)

Return the numeric precision value stored on the graph. This is used
as a default metadata value and as a parameter for certain distance
computations (e.g. `shape_distance`).
"""
function precision(g::DynamicGeometricGraph)
    return g.precision
end

"""
get_coords(g)

Return a tuple `(P, vs)` where `P` is a matrix with vertex coordinates as
columns (sorted by vertex id) and `vs` is the vector of vertex ids in the
same order. Useful for plotting or exporting coordinates in bulk.
"""
function get_coords(g::DynamicGeometricGraph)
    vs = sort(vertices(g)|>collect)
    Ps=[get_vertex_coords(g, v) for v in vs]
    return permutedims(hcat(Ps...)), vs
end


"""
geomfilter(wm, wM)

Internal helper to combine a minimum (`wm`) and maximum (`wM`) weight into a
single filtered value. Used when deriving an aggregate edge weight. If
`wm == 0` the function returns `0`.
"""
function geomfilter(wm, wM)
    if wm > wM
        t = wm
        wm = wM
        wM = t
    end
    if wm == 0
        return 0
    end
    return (sqrt(wm*wM)) * (wm / wM )
end

"""
weighted_degree(g, vertex)

Return the sum of weights for edges incident to `vertex`. If the vertex has
no incident edges, `0` is returned and a warning is emitted.
"""
function weighted_degree(g, vertex)
    ies = incident_edges(g, vertex)|>collect
    if length(ies) == 0
        @debug "0-degree for $vertex"
        return 0
    end
    return sum(w for (_, _, w) in ies)
end

"""
get_neighbour_coords(g, vertex; sort=false)

Return a pair `(coords, ids)` with coordinates and vertex ids of neighbours
of `vertex`. If `sort=true` the neighbour coordinates are ordered
clockwise around the query vertex (uses `sort_clockwise_indices`).
"""
function get_neighbour_coords(g::DynamicGeometricGraph, vertex::Int, sort=false)
    ies = incident_edges(g, vertex)
    vs = [u for (_, u, _) in ies ]
    crds = [get_vertex_coords(g, u) for u in vs]
    if sort
        idx = sort_clockwise_indices(crds)
        return crds[idx], vs[idx]
    else
        crds, vs
    end
end

"""
sort_clockwise_indices(points)

Return an index permutation that orders 2D `points` clockwise around their
centroid. Intended for small neighbour lists; uses `atan` for angle sort.
"""
function sort_clockwise_indices(points::Vector{SVector{2, Float64}})
    # Calculate centroid
    centroid = sum(points) / length(points)
    
    # Get the permutation of indices that would sort the points clockwise
    indices = sortperm(points, by=p -> atan(p[2] - centroid[2], p[1] - centroid[1]), rev=true)
    return indices
end


"""
    graph_edits(SpatialGraph)

    Return the number of graph edits done on this graph since its construction or copy.
    Adding, removing a vertex count as one. 
    Adding or removing an edge count as one.
    Removing a vertex with k adjacent edges is 1+k edits.
"""
function graph_edits(sg)
    return sg.counter[2]
end


Base.copy(g::DynamicGeometricGraph)=DynamicGeometricGraphs.freeze(g)

"""
update_graph(g)

Internal helper to record that the graph has been edited. Increments the
internal edit counter. Users typically do not call this directly; it is
invoked by mutating operations such as `add_vertex!` or `add_edge!`.
"""
function update_graph(g)
    g.counter[2] += 1
end



"""
euclid(a, b)

Euclidean distance between two points. Several method overloads are
provided: scalar inputs, 2D `SVector{2,T}` and 3D `SVector{3,T}`. Returns a
non-negative real value.
"""
function euclid(a::SVector{2, T}, b::SVector{2,T}) where {T <: Real}
    sqrt(sum((a .- b).^2))
end

function euclid(a::T, b::T) where {T <: Real}
    sqrt(sum((a - b)^2))
end


function euclid(a::SVector{3, T}, b::SVector{3,T}) where {T <: Real}
    sqrt(sum((a .- b).^2))
end


# ============================================================================
# Graph transformation and generation functions
# ============================================================================

"""
    update_coord!(g::DynamicGeometricGraph{2, T}, old::SVector{2, T}, new::SVector{2, T}) where T<:Real

Update vertex coordinates from old to new position. Modifies the graph in place.

# Arguments
- `g`: The graph to modify
- `old`: The current coordinates of the vertex
- `new`: The new coordinates for the vertex

# Returns
- The modified graph `g`
"""
function update_coord!(g::DynamicGeometricGraph{2, T}, old::SVector{2, T}, new::SVector{2, T}) where T<:Real
    if old == new
        @debug "Old coordinates are new, not doing anything"
        return g
    end
    # Find the vertex the coordinates correspond to
    vid = g.lookupvertices[old]
    # Check that the new coordinates don't by accident already are used
    if new in keys(g.lookupvertices)
        @error "Invalid new coordinates $(g.lookupvertices[new])"
        # In this case, we have two options, we assume that the caller knows what they're doing
        # But then we need to 'fuse' vertices, and will break planarity
        # This is costly, we need to drop the vertex, and reassing all edges
        # The caller instead should be warned, and the caller can check by
        # lookup(new)
        return g
    end
    # The normal case, update the coordinates
    g.lookupvertices[new] = vid
    # index -> coord
    g.vertices[vid] = new
    # Don't forget to remove stale coordinates
    delete!(g.lookupvertices, old)
    update_graph(g)
    return g
end


"""
    find_nearest(g::DynamicGeometricGraph{N, T, F}, v::SVector{N, T}) where {N, T, F}

Find the nearest vertex in the graph to the given coordinates.

# Arguments
- `g`: The graph to search
- `v`: The coordinates to find the nearest vertex to

# Returns
- The vertex ID of the nearest vertex
"""
function find_nearest(g::DynamicGeometricGraph{N, T, F}, v::SVector{N, T}) where {N, T, F}
    if haskey(g.lookupvertices, v)
        @debug "Found an exact match for $v"
        return g.lookupvertices[v]
    else
        @debug "Finding nearest vertex for $v"
        isempty(g.vertices) && error("Graph has no vertices")
        best_id = -1
        best_d = typemax(T)
        for (vi, vc) in g.vertices
            d = g.distance(v, vc)
            if d < best_d
                best_d = d
                best_id = vi
            end
        end
        return best_id
    end
end


"""
    generate_hexagonal_points(j::Float64)

Generate 7 points in a hexagonal arrangement: one center point at origin
and 6 points arranged in a hexagon around it.

# Arguments
- `j`: The distance from the center to each outer point

# Returns
- A vector of 7 tuples representing (x, y) coordinates
"""
function generate_hexagonal_points(j::Float64)
    # The center point is at the origin
    center = (0.0, 0.0)

    # Initialize a vector to hold the coordinates, starting with the center
    locations = [center]

    # Calculate the coordinates of the 6 outer points
    for i in 0:5
        angle = i * pi / 3
        x = j * cos(angle)
        y = j * sin(angle)
        push!(locations, (x, y))
    end

    return locations
end


"""
    generate_hub_spoke_graph(k, n::Int)

Generate a hub-spoke graph with a central hub and n spokes radiating outward.

# Arguments
- `k`: The distance from hub to spoke endpoints
- `n`: The number of spokes (must be between 1 and 6)

# Returns
- A `DynamicGeometricGraph` with the hub-spoke structure
"""
function generate_hub_spoke_graph(k, n::Int)
    g=DynamicGeometricGraph{2, Float64}()
    # Validate that the number of spokes is within the specified range (1-6).
    if !(1 <= n <= 6)
        error("The number of connected vertices (n) must be between 1 and 6.")
    end

    # Create the hub at the origin and add it to the graph.
    hub_coords = SVector{2, Float64}(0.0, 0.0)
    hub_index = DynamicGeometricGraphs.add_vertex!(g, hub_coords)

    # Define the angular step between spokes.
    angle_step = pi / 3

    # Iterate to create each spoke.
    for i in 0:(n-1)
        # Calculate the angle for the current spoke. A negative sign ensures clockwise placement.
        angle = -i * angle_step

        # Determine the coordinates of the spoke's endpoint using polar to Cartesian conversion.
        x = k * cos(angle)
        y = k * sin(angle)
        spoke_coords = SVector{2, Float64}(x, y)

        # Add the spoke vertex to the graph.
        spoke_index = DynamicGeometricGraphs.add_vertex!(g, spoke_coords)

        # Connect the new spoke vertex to the central hub.
        DynamicGeometricGraphs.add_edge!(g, hub_index, spoke_index)
    end

    return g
end


"""
    translate_graph(g::DynamicGeometricGraph{2, T}, translation_vector::SVector{2, T}) where T<:Real

Creates a translated copy of the graph by adding a translation vector to all vertices.

# Arguments
- `g`: The source graph to translate
- `translation_vector`: The vector to add to all vertex coordinates

# Returns
- A new `DynamicGeometricGraph` object with all vertices translated
"""
function translate_graph(g::DynamicGeometricGraph{2, T}, translation_vector::SVector{2, T}) where T<:Real
    g_new = freeze(g)

    crds, verts = get_coords(g_new)

    for (i, v_id) in enumerate(verts)
        old_coords = SVector{2, T}(crds[i, 1], crds[i, 2])
        new_coords = old_coords + translation_vector
        update_coord!(g_new, old_coords, new_coords)
    end

    return g_new
end


"""
    scale_graph(g::DynamicGeometricGraph{2, T}, factor::T) where T<:Real

Creates a scaled copy of the graph relative to the origin (0,0).

# Arguments
- `g`: The source graph to scale
- `factor`: The scaling factor to apply to all coordinates

# Returns
- A new `DynamicGeometricGraph` object with all vertices scaled
"""
function scale_graph(g::DynamicGeometricGraph{2, T}, factor::T) where T<:Real
    g_new = freeze(g)

    crds, verts = get_coords(g_new)

    for (i, v_id) in enumerate(verts)
        old_coords = SVector{2, T}(crds[i, 1], crds[i, 2])
        new_coords = old_coords * factor
        update_coord!(g_new, old_coords, new_coords)
    end

    return g_new
end


"""
    rotate_graph(g::DynamicGeometricGraph{2, T}, angle::T) where T<:Real

Creates a rotated copy of the graph around the origin (0,0).

# Arguments
- `g`: The source graph to rotate
- `angle`: The rotation angle in radians

# Returns
- A new `DynamicGeometricGraph` object with all vertices rotated
"""
function rotate_graph(g::DynamicGeometricGraph{2, T}, angle::T) where T<:Real
    g_new = freeze(g)
    c, s = cos(angle), sin(angle)
    rotation_matrix = @SMatrix [c -s; s c]

    crds, verts = get_coords(g_new)

    for (i, v_id) in enumerate(verts)
        old_coords = SVector{2, T}(crds[i, 1], crds[i, 2])
        new_coords = rotation_matrix * old_coords
        update_coord!(g_new, old_coords, new_coords)
    end

    return g_new
end


"""
    merge_graphs(f, g)

Merge two graphs by copying all vertices and edges from g into a copy of f.

# Arguments
- `f`: The first graph (will be copied)
- `g`: The second graph (vertices and edges will be added to copy of f)

# Returns
- A new graph containing all vertices and edges from both graphs
"""
function merge_graphs(f, g)
    h = freeze(f)
    b = Bijections.Bijection{Int, Int}()
    for v in vertices(g)
        c = get_vertex_coords(g, v)
        m = get_metadata(g, v)
        vi=DynamicGeometricGraphs.add_vertex!(h, c; meta=m)
        # set_metadata!(h, vi, m)
        b[v]=vi
    end
    for e in edges(g)
        u = src(e)
        v = dst(e)
        ui = b[u]
        vi = b[v]
        add_edge!(h, ui, vi)
    end
    h
end


"""
    refgraph(ctrspoke=300.0, motifdistance=100; pattern=[1,2,3,4,5,6])

Generate a reference graph with a central hub-spoke structure surrounded by
hexagonally arranged motif graphs.

# Arguments
- `ctrspoke`: Distance from center to the hexagonal arrangement points (default: 300.0)
- `motifdistance`: Distance parameter for the hub-spoke motifs (default: 100)
- `pattern`: Array specifying the number of spokes for each of the 6 motif positions (default: [1,2,3,4,5,6])

# Returns
- A `DynamicGeometricGraph` containing the complete reference graph structure
"""
function refgraph(ctrspoke=300.0, motifdistance=100; pattern=[1,2,3,4,5,6])
    pts = generate_hexagonal_points(ctrspoke)
    vecs = [SVector{2, Float64}.((pts[i] .- pts[1])...) for i in 2:length(pts)]
    gb = DynamicGeometricGraph{2, Float64}()
    ctr = SVector{2, Float64}(0, 0)
    id=add_vertex!(gb, ctr)
    graphs = [generate_hub_spoke_graph(motifdistance, n) for n in pattern]
    for (i,g) in enumerate(graphs)
        tv = translate_graph(g, vecs[i])
        crds, vs  = get_coords(tv)
        distances = [LinearAlgebra.norm(crds[i,:] - ctr) for i in 1:size(crds, 1)]
        _ , row_id = findmin(distances)
        sv = crds[row_id, :]
        gb = merge_graphs(gb, translate_graph(g, vecs[i]))
        v = find_nearest(gb, sv)
        add_edge!(gb, 1, v)
    end
    return gb
end


end # module
