# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Ben Cardoen <b.cardoen@bham.ac.uk>

"""
    generate_wavefront_viz.jl

Generate an animated GIF showing a DNA-like double helix with dynamic cross-edges.
Two spiral strands have permanent backbone edges (degree 2). Cross-edges ("rungs")
between the strands appear and disappear over time.

Usage:
    julia --project=. scripts/generate_wavefront_viz.jl
"""

using DynamicGeometricGraphs
using Graphs
using StaticArrays
using CairoMakie
using Random

Random.seed!(42)

# ---------------------------------------------------------------------------
# 1. Build double-helix points
# ---------------------------------------------------------------------------
n_per_strand = 40                # vertices per strand (80 total)
helix_r      = 60.0              # helix radius
helix_h      = 300.0             # total height
n_turns      = 3.0               # number of full turns

strand_a = Vector{SVector{3,Float64}}(undef, n_per_strand)
strand_b = Vector{SVector{3,Float64}}(undef, n_per_strand)

for i in 1:n_per_strand
    t     = (i - 1) / (n_per_strand - 1)   # 0..1
    angle = 2π * n_turns * t
    z     = helix_h * t
    strand_a[i] = SVector(helix_r * cos(angle), helix_r * sin(angle), z)
    strand_b[i] = SVector(helix_r * cos(angle + π), helix_r * sin(angle + π), z)
end

# ---------------------------------------------------------------------------
# 2. Build the graph: static vertices + backbone edges
# ---------------------------------------------------------------------------
g = DynamicGeometricGraph{3,Float64}()

ids_a = Int[]
ids_b = Int[]

for i in 1:n_per_strand
    push!(ids_a, add_vertex!(g, strand_a[i]))
    push!(ids_b, add_vertex!(g, strand_b[i]))
end

# Backbone edges (permanent, degree 2 chain per strand)
for i in 1:(n_per_strand - 1)
    add_edge!(g, ids_a[i], ids_a[i+1])
    add_edge!(g, ids_b[i], ids_b[i+1])
end

# ---------------------------------------------------------------------------
# 3. Cross-edge candidates: pairs at same height (rungs)
# ---------------------------------------------------------------------------
rung_pairs = [(ids_a[i], ids_b[i]) for i in 1:n_per_strand]

# ---------------------------------------------------------------------------
# 4. Animate: sweep a window of active rungs up and down
# ---------------------------------------------------------------------------
n_frames     = 40
rung_window  = 8       # how many rungs are active at once
sweep_speed  = 1.5     # rungs per frame

struct FrameData
    # strand coords (static)
    ax::Vector{Float64}; ay::Vector{Float64}; az::Vector{Float64}
    bx::Vector{Float64}; by::Vector{Float64}; bz::Vector{Float64}
    # active rung endpoints
    rx1::Vector{Float64}; ry1::Vector{Float64}; rz1::Vector{Float64}
    rx2::Vector{Float64}; ry2::Vector{Float64}; rz2::Vector{Float64}
    title::String
end

frames = FrameData[]

for fr in 1:n_frames
    # Remove old rungs
    for (u, v) in rung_pairs
        if Graphs.has_edge(g, u, v)
            rem_edge!(g, u, v)
        end
    end

    # Compute which rungs are active: a window that sweeps up
    center = 1.0 + (fr - 1) * sweep_speed
    active_rungs = Int[]
    for i in 1:n_per_strand
        if abs(i - center) < rung_window / 2
            push!(active_rungs, i)
        end
    end

    # Add active rungs
    for i in active_rungs
        u, v = rung_pairs[i]
        add_edge!(g, u, v)
    end

    # Snapshot
    sa_x = [strand_a[i][1] for i in 1:n_per_strand]
    sa_y = [strand_a[i][2] for i in 1:n_per_strand]
    sa_z = [strand_a[i][3] for i in 1:n_per_strand]
    sb_x = [strand_b[i][1] for i in 1:n_per_strand]
    sb_y = [strand_b[i][2] for i in 1:n_per_strand]
    sb_z = [strand_b[i][3] for i in 1:n_per_strand]

    rx1, ry1, rz1 = Float64[], Float64[], Float64[]
    rx2, ry2, rz2 = Float64[], Float64[], Float64[]
    for i in active_rungs
        a, b = strand_a[i], strand_b[i]
        push!(rx1, a[1]); push!(ry1, a[2]); push!(rz1, a[3])
        push!(rx2, b[1]); push!(ry2, b[2]); push!(rz2, b[3])
    end

    title = "nv=$(nv(g))  ne=$(ne(g))  rungs=$(length(active_rungs))"
    push!(frames, FrameData(sa_x, sa_y, sa_z, sb_x, sb_y, sb_z,
                            rx1, ry1, rz1, rx2, ry2, rz2, title))
end

# ---------------------------------------------------------------------------
# 5. Render
# ---------------------------------------------------------------------------
println("Rendering $(length(frames)) frames ...")

output_path = joinpath(@__DIR__, "..", "docs", "dynamic_graph_operations.gif")

fig = Figure(size=(420, 520), backgroundcolor=:white)

pad = 20.0
all_x = vcat([p[1] for p in strand_a], [p[1] for p in strand_b])
all_y = vcat([p[2] for p in strand_a], [p[2] for p in strand_b])
all_z = vcat([p[3] for p in strand_a], [p[3] for p in strand_b])

record(fig, output_path, 1:length(frames); framerate=6) do i
    fd = frames[i]
    empty!(fig)

    ax = Axis3(fig[1, 1],
               title = fd.title,
               titlesize = 14,
               aspect = (1, 1, 2.0),
               azimuth = 0.3π + 0.04π * i,   # slow rotation
               elevation = 0.15π,
               backgroundcolor = :grey98)
    xlims!(ax, minimum(all_x) - pad, maximum(all_x) + pad)
    ylims!(ax, minimum(all_y) - pad, maximum(all_y) + pad)
    zlims!(ax, minimum(all_z) - pad, maximum(all_z) + pad)

    # Backbone strands
    lines!(ax, fd.ax, fd.ay, fd.az; color=:steelblue, linewidth=2.0)
    lines!(ax, fd.bx, fd.by, fd.bz; color=:indianred, linewidth=2.0)

    # Vertices on strands
    scatter!(ax, fd.ax, fd.ay, fd.az;
             color=:steelblue, markersize=6, strokewidth=0.5, strokecolor=:white)
    scatter!(ax, fd.bx, fd.by, fd.bz;
             color=:indianred, markersize=6, strokewidth=0.5, strokecolor=:white)

    # Cross-rungs
    for k in eachindex(fd.rx1)
        lines!(ax,
               [fd.rx1[k], fd.rx2[k]],
               [fd.ry1[k], fd.ry2[k]],
               [fd.rz1[k], fd.rz2[k]];
               color=(:goldenrod, 0.8), linewidth=2.5)
    end
end

println("Saved to $output_path")
