# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Ben Cardoen <b.cardoen@bham.ac.uk>

"""
    generate_metric_comparison_viz.jl

Demonstrate different distance metrics on the same graph structure.
"""

using DynamicGeometricGraphs
using Graphs
using StaticArrays
using GLMakie
using LinearAlgebra

println("Generating metric comparison...")

# Define kernel distance
function kernel_distance(a, b; sigma=1.0)
    sq_dist = sum((a .- b).^2)
    return 1.0 - exp(-sq_dist / (2 * sigma^2))
end

# Build graphs
g_euclidean = DynamicGeometricGraph{3, Float64}()
g_kernel = DynamicGeometricGraph{3, Float64}(distfun=(a,b)->kernel_distance(a,b, sigma=1.5))

coords = [SVector{3, Float64}(0.0, 0.0, 0.0),
          SVector{3, Float64}(2.0, 0.0, 0.0),
          SVector{3, Float64}(1.0, 2.0, 0.0),
          SVector{3, Float64}(1.0, 1.0, 2.0)]

v_euclidean = [add_vertex!(g_euclidean, c) for c in coords]
v_kernel = [add_vertex!(g_kernel, c) for c in coords]

for (i, j) in [(1,2), (2,3), (3,1), (1,4), (2,4), (3,4)]
    add_edge!(g_euclidean, v_euclidean[i], v_euclidean[j])
    add_edge!(g_kernel, v_kernel[i], v_kernel[j])
end

# Get weights
function get_edge_weights(g)
    weights = Dict{Tuple{Int,Int}, Float64}()
    for (u, neighbors) in g.edges
        for v in neighbors
            if u < v
                weights[(u, v)] = edge_weight(g, u, v)
            end
        end
    end
    return weights
end

weights_euc = get_edge_weights(g_euclidean)
weights_kern = get_edge_weights(g_kernel)

println("Euclidean weights: ", sort(collect(weights_euc)))
println("Kernel weights: ", sort(collect(weights_kern)))

wmin = minimum(vcat(collect(values(weights_euc)), collect(values(weights_kern))))
wmax = maximum(vcat(collect(values(weights_euc)), collect(values(weights_kern))))

# Create static side-by-side comparison
fig = Figure(size=(1200, 550))

Label(fig[0,:], "Same Geometry, Different Metrics", fontsize=26, font=:bold)

ax1 = Axis3(fig[1,1], title="Euclidean Distance", azimuth=0.3π, elevation=0.15π)
ax2 = Axis3(fig[1,2], title="Kernel Distance (σ=1.5)", azimuth=0.3π, elevation=0.15π)

M, _ = get_coords(g_euclidean)
x_c = [M[i,1] for i in 1:size(M,1)]
y_c = [M[i,2] for i in 1:size(M,1)]
z_c = [M[i,3] for i in 1:size(M,1)]

# Euclidean (left)
for ((u,v), w) in weights_euc
    u_c, v_c = g_euclidean.vertices[u], g_euclidean.vertices[v]
    lines!(ax1, [u_c[1], v_c[1]], [u_c[2], v_c[2]], [u_c[3], v_c[3]],
           color=(w-wmin)/(wmax-wmin), colormap=:viridis, colorrange=(0,1), linewidth=6)
    mid = (u_c .+ v_c) ./ 2
    text!(ax1, "$(round(w, digits=2))", position=Point3f(mid...), fontsize=14)
end
scatter!(ax1, x_c, y_c, z_c, color=:white, markersize=22,
        strokewidth=2, strokecolor=:steelblue)

# Kernel (right)
for ((u,v), w) in weights_kern
    u_c, v_c = g_kernel.vertices[u], g_kernel.vertices[v]
    lines!(ax2, [u_c[1], v_c[1]], [u_c[2], v_c[2]], [u_c[3], v_c[3]],
           color=(w-wmin)/(wmax-wmin), colormap=:viridis, colorrange=(0,1), linewidth=6)
    mid = (u_c .+ v_c) ./ 2
    text!(ax2, "$(round(w, digits=2))", position=Point3f(mid...), fontsize=14)
end
scatter!(ax2, x_c, y_c, z_c, color=:white, markersize=22,
        strokewidth=2, strokecolor=:darkorange)

xlims!(ax1, -0.5, 2.5); ylims!(ax1, -0.5, 2.5); zlims!(ax1, -0.5, 2.5)
xlims!(ax2, -0.5, 2.5); ylims!(ax2, -0.5, 2.5); zlims!(ax2, -0.5, 2.5)

Colorbar(fig[1,3], limits=(wmin, wmax), colormap=:viridis, label="Weight")

save("docs/metric_comparison.png", fig)
println("✓ Saved to docs/metric_comparison.png")
