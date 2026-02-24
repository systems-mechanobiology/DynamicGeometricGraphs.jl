# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Ben Cardoen <b.cardoen@bham.ac.uk>

"""
    generate_dynamic_viz.jl

Generate an animated GIF demonstrating dynamic vertex and edge operations in DynamicGeometricGraphs.
Shows add/remove operations for both vertices and edges in 3D.

Usage:
    julia --project=. scripts/generate_dynamic_viz.jl
"""

using DynamicGeometricGraphs
using Graphs
using StaticArrays
using GLMakie

println("Generating 3D dynamic graph visualization...")

# Create frames for the animation
frames = []

# Frame 1: Empty graph
g = DynamicGeometricGraph{3, Float64}()
push!(frames, (freeze(g), "Step 1: Empty graph"))

# Frame 2-4: Build a tetrahedron base (4 vertices)
v1 = add_vertex!(g, SVector{3, Float64}(0.0, 0.0, 0.0))
push!(frames, (freeze(g), "Step 2: Add vertex v1 (0,0,0)"))

v2 = add_vertex!(g, SVector{3, Float64}(2.0, 0.0, 0.0))
push!(frames, (freeze(g), "Step 3: Add vertex v2 (2,0,0)"))

v3 = add_vertex!(g, SVector{3, Float64}(1.0, 2.0, 0.0))
push!(frames, (freeze(g), "Step 4: Add vertex v3 (1,2,0)"))

# Frame 5: Create base triangle edges
add_edge!(g, v1, v2)
add_edge!(g, v2, v3)
add_edge!(g, v3, v1)
push!(frames, (freeze(g), "Step 5: Form base triangle"))

# Frame 6: Add apex vertex
v4 = add_vertex!(g, SVector{3, Float64}(1.0, 1.0, 2.0))
push!(frames, (freeze(g), "Step 6: Add apex v4 (1,1,2)"))

# Frame 7-9: Connect apex to form tetrahedron
add_edge!(g, v1, v4)
push!(frames, (freeze(g), "Step 7: Connect v1-v4"))

add_edge!(g, v2, v4)
push!(frames, (freeze(g), "Step 8: Connect v2-v4"))

add_edge!(g, v3, v4)
push!(frames, (freeze(g), "Step 9: Complete tetrahedron"))

# Frame 10: Add a new vertex outside
v5 = add_vertex!(g, SVector{3, Float64}(3.5, 1.0, 1.0))
push!(frames, (freeze(g), "Step 10: Add external vertex v5"))

# Frame 11: Connect to multiple vertices
add_edge!(g, v2, v5)
add_edge!(g, v4, v5)
push!(frames, (freeze(g), "Step 11: Connect v5 to structure"))

# Frame 12: Add another external vertex
v6 = add_vertex!(g, SVector{3, Float64}(-1.0, 1.0, 1.0))
add_edge!(g, v1, v6)
add_edge!(g, v3, v6)
push!(frames, (freeze(g), "Step 12: Add v6 and connections"))

# Frame 13: Remove an edge
rem_edge!(g, v1, v2)
push!(frames, (freeze(g), "Step 13: Remove edge v1-v2"))

# Frame 14: Remove a vertex (and its edges)
rem_vertex!(g, v5)
push!(frames, (freeze(g), "Step 14: Remove vertex v5"))

# Frame 15: Remove another edge
rem_edge!(g, v3, v4)
push!(frames, (freeze(g), "Step 15: Remove edge v3-v4"))

# Frame 16-17: Hold final state
push!(frames, (freeze(g), "Step 16: Final graph structure"))
push!(frames, (freeze(g), "Step 17: Complete"))

println("Rendering $(length(frames)) frames in 3D...")

# Record the animation
output_path = "docs/dynamic_graph_operations.gif"
fig = Figure(size=(600, 600))

record(fig, output_path, 1:length(frames), framerate=2) do i
    graph, title = frames[i]
    println("  Rendering frame $i/$(length(frames)): $title")
    
    empty!(fig)
    ax = Axis3(fig[1, 1], 
               title=title,
               xlabel="x", ylabel="y", zlabel="z",
               aspect=(1, 1, 1),
               azimuth=0.3π,
               elevation=0.15π)
    
    # Get coordinates
    if nv(graph) > 0
        M, vids = get_coords(graph)
        
        # Plot edges first (behind vertices)
        for (u_vid, neighbors) in graph.edges
            if haskey(graph.vertices, u_vid)
                u_coords = graph.vertices[u_vid]
                for v_vid in neighbors
                    if haskey(graph.vertices, v_vid) && u_vid < v_vid
                        v_coords = graph.vertices[v_vid]
                        lines!(ax, [u_coords[1], v_coords[1]], 
                                  [u_coords[2], v_coords[2]],
                                  [u_coords[3], v_coords[3]],
                                  color=(:gray, 0.6), linewidth=3)
                    end
                end
            end
        end
        
        # Plot vertices on top - convert to regular arrays
        x_coords = [M[i, 1] for i in 1:size(M, 1)]
        y_coords = [M[i, 2] for i in 1:size(M, 1)]
        z_coords = [M[i, 3] for i in 1:size(M, 1)]
        scatter!(ax, x_coords, y_coords, z_coords, 
                color=:steelblue, markersize=20,
                strokewidth=2, strokecolor=:white)
    end
    
    # Set consistent limits
    xlims!(ax, -1.5, 4.0)
    ylims!(ax, -0.5, 2.5)
    zlims!(ax, -0.5, 2.5)
end

println("✓ Animation saved to: $output_path")
println("\nDemonstrates:")
println("  • 3D geometric graph with N-dimensional coordinates")
println("  • Building complex 3D structures (tetrahedron + extensions)")
println("  • Adding vertices in 3D space")
println("  • Creating edges to form polyhedra")
println("  • Extending structures with new vertices")
println("  • Removing edges dynamically")
println("  • Removing vertices dynamically")
println("  • Graph remains consistent throughout all operations")
println("  • $(length(frames)) steps showing progressive construction and modification")
