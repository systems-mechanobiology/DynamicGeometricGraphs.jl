# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Ben Cardoen <b.cardoen@bham.ac.uk>

using Test
using Graphs
using StaticArrays
using DynamicGeometricGraphs
import Base.CoreLogging: NullLogger, with_logger

@testset "DynamicGeometricGraph" begin
    @testset "Graphs.jl API" begin
        @testset "Basic construction and properties" begin
            g = DynamicGeometricGraph{2, Float64}()
            @test nv(g) == 0
            @test ne(g) == 0
            @test isempty(collect(vertices(g)))
            @test isempty(collect(edges(g)))
        end

        @testset "Vertex operations" begin
            g = DynamicGeometricGraph{2, Float64}()

            # Add vertices
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(1.0, 0.0))
            v3 = add_vertex!(g, SVector{2, Float64}(1.0, 1.0))

            @test nv(g) == 3
            @test v1 == 1
            @test v2 == 2
            @test v3 == 3

            # Test has_vertex (via Graphs API)
            @test Graphs.has_vertex(g, v1)
            @test Graphs.has_vertex(g, v2)
            @test Graphs.has_vertex(g, v3)
            @test !Graphs.has_vertex(g, 999)

            # Test vertices iterator
            verts = collect(vertices(g))
            @test length(verts) == 3
            @test v1 ∈ verts
            @test v2 ∈ verts
            @test v3 ∈ verts

            # Add duplicate vertex (should return existing index)
            v1_dup = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            @test v1_dup == v1
            @test nv(g) == 3
        end

        @testset "Edge operations" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(1.0, 0.0))
            v3 = add_vertex!(g, SVector{2, Float64}(1.0, 1.0))

            # Initially no edges
            @test ne(g) == 0
            @test !Graphs.has_edge(g, v1, v2)

            # Add edges
            add_edge!(g, v1, v2)
            @test ne(g) == 1
            @test Graphs.has_edge(g, v1, v2)
            @test Graphs.has_edge(g, v2, v1)  # Undirected

            add_edge!(g, v2, v3)
            add_edge!(g, v1, v3)
            @test ne(g) == 3

            # Test edges iterator
            edge_list = collect(edges(g))
            @test length(edge_list) == 3
        end

        @testset "Neighbors" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(1.0, 0.0))
            v3 = add_vertex!(g, SVector{2, Float64}(1.0, 1.0))

            add_edge!(g, v1, v2)
            add_edge!(g, v1, v3)

            # Test outneighbors
            out_nbrs = collect(outneighbors(g, v1))
            @test length(out_nbrs) == 2
            @test v2 ∈ out_nbrs
            @test v3 ∈ out_nbrs

            # Test inneighbors (same as outneighbors for undirected)
            in_nbrs = collect(inneighbors(g, v1))
            @test length(in_nbrs) == 2
            @test v2 ∈ in_nbrs
            @test v3 ∈ in_nbrs
        end

        @testset "Vertex removal" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(1.0, 0.0))
            v3 = add_vertex!(g, SVector{2, Float64}(1.0, 1.0))

            add_edge!(g, v1, v2)
            add_edge!(g, v1, v3)
            add_edge!(g, v2, v3)

            @test nv(g) == 3
            @test ne(g) == 3

            # Remove vertex
            rem_vertex!(g, v1)
            @test nv(g) == 2
            @test !Graphs.has_vertex(g, v1)
            @test ne(g) == 1  # Only v2-v3 edge remains
            @test Graphs.has_edge(g, v2, v3)
            @test !Graphs.has_edge(g, v1, v2)
            @test !Graphs.has_edge(g, v1, v3)
        end

        @testset "Edge removal (rem_edge!)" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(1.0, 0.0))
            v3 = add_vertex!(g, SVector{2, Float64}(1.0, 1.0))
            v4 = add_vertex!(g, SVector{2, Float64}(0.0, 1.0))

            # Add edges to form a square
            add_edge!(g, v1, v2)
            add_edge!(g, v2, v3)
            add_edge!(g, v3, v4)
            add_edge!(g, v4, v1)

            @test ne(g) == 4
            @test Graphs.has_edge(g, v1, v2)
            @test Graphs.has_edge(g, v2, v1)  # Undirected

            # Remove one edge
            rem_edge!(g, v1, v2)
            @test ne(g) == 3
            @test !Graphs.has_edge(g, v1, v2)
            @test !Graphs.has_edge(g, v2, v1)  # Both directions removed

            # Other edges should still exist
            @test Graphs.has_edge(g, v2, v3)
            @test Graphs.has_edge(g, v3, v4)
            @test Graphs.has_edge(g, v4, v1)

            # Remove another edge
            rem_edge!(g, v3, v4)
            @test ne(g) == 2
            @test !Graphs.has_edge(g, v3, v4)

            # Remove same edge again (should be no-op)
            rem_edge!(g, v3, v4)
            @test ne(g) == 2  # Count unchanged

            # Test removing edge with non-existent vertex (should warn, not error)
            rem_edge!(g, v1, 999)
            @test ne(g) == 2  # Graph unchanged

            # Remove all remaining edges
            rem_edge!(g, v2, v3)
            rem_edge!(g, v4, v1)
            @test ne(g) == 0
            @test nv(g) == 4  # Vertices still exist
        end

        @testset "3D graph" begin
            g = DynamicGeometricGraph{3, Float64}()
            v1 = add_vertex!(g, SVector{3, Float64}(0.0, 0.0, 0.0))
            v2 = add_vertex!(g, SVector{3, Float64}(1.0, 0.0, 0.0))
            v3 = add_vertex!(g, SVector{3, Float64}(0.0, 1.0, 0.0))

            @test nv(g) == 3
            add_edge!(g, v1, v2)
            add_edge!(g, v2, v3)
            @test ne(g) == 2
            @test Graphs.has_vertex(g, v1)
            @test Graphs.has_edge(g, v1, v2)
        end
    end

    @testset "Package API" begin
        @testset "Coordinate operations" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(1.0, 2.0))
            v2 = add_vertex!(g, SVector{2, Float64}(3.0, 4.0))
            v3 = add_vertex!(g, SVector{2, Float64}(5.0, 6.0))

            # Test get_vertex_coords
            coords = get_vertex_coords(g, v1)
            @test coords == SVector{2, Float64}(1.0, 2.0)

            # Test get_vertex_idx
            idx = get_vertex_idx(g, SVector{2, Float64}(3.0, 4.0))
            @test idx == v2

            # Test has_vertex with SVector (package function)
            @test DynamicGeometricGraphs.has_vertex(g, SVector{2, Float64}(1.0, 2.0))
            @test !DynamicGeometricGraphs.has_vertex(g, SVector{2, Float64}(99.0, 99.0))

            # Test get_coords
            coord_matrix, vertex_ids = get_coords(g)
            @test size(coord_matrix) == (3, 2)
            @test length(vertex_ids) == 3
        end

        @testset "Metadata operations" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))

            # Initial metadata should be precision-based
            meta = get_metadata(g, v1)
            @test meta == SVector{2, Float64}(1.0, 1.0)

            # Set new metadata
            new_meta = SVector{2, Float64}(5.0, 10.0)
            set_metadata!(g, v1, new_meta)
            @test get_metadata(g, v1) == new_meta

            # Test non-existent vertex
            @test get_metadata(g, 999) === nothing
        end

        @testset "Edge weight and distance" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(3.0, 4.0))
            add_edge!(g, v1, v2)

            # Test spatial edge weight
            w = edge_weight(g, v1, v2; adjacency=:spatial)
            @test w ≈ 5.0  # 3-4-5 triangle

            # Test with custom distance function
            custom_dist(a, b) = sum(abs.(a .- b))  # Manhattan distance
            w_custom = edge_weight(g, v1, v2; adjacency=:spatial, distancefun=custom_dist)
            @test w_custom == 7.0

            # Test shape-based edge weight (experimental feature, suppress warning)
            set_metadata!(g, v1, SVector{2, Float64}(1.0, 1.0))
            set_metadata!(g, v2, SVector{2, Float64}(2.0, 2.0))
            with_logger(NullLogger()) do
                w_shape = edge_weight(g, v1, v2; adjacency=:shape)
                @test 0.0 <= w_shape <= 2.0
            end
        end

        @testset "Euclidean distance" begin
            # 2D
            a2 = SVector{2, Float64}(0.0, 0.0)
            b2 = SVector{2, Float64}(3.0, 4.0)
            @test euclid(a2, b2) ≈ 5.0

            # 3D
            a3 = SVector{3, Float64}(0.0, 0.0, 0.0)
            b3 = SVector{3, Float64}(1.0, 2.0, 2.0)
            @test euclid(a3, b3) ≈ 3.0

            # Scalar
            @test euclid(5.0, 8.0) ≈ 3.0
        end

        @testset "Shape distance and similarity" begin
            # Similarity
            @test similarity(1.0, 1.0) ≈ 1.0  # Identical
            @test 0.0 < similarity(1.0, 2.0) < 1.0  # Different

            # Shape distance
            d = shape_distance(1.0, 1.0)
            @test d ≈ 0.0  # Identical
            d2 = shape_distance(1.0, 5.0)
            @test d2 > 0.0  # Different
        end

        @testset "Incident edges and neighbors" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(1.0, 0.0))
            v3 = add_vertex!(g, SVector{2, Float64}(1.0, 1.0))
            add_edge!(g, v1, v2)
            add_edge!(g, v1, v3)

            # Test incident_edges
            ies = collect(incident_edges(g, v1))
            @test length(ies) == 2

            # Test get_neighbour_coords without sorting
            coords, verts = get_neighbour_coords(g, v1, false)
            @test length(coords) == 2
            @test length(verts) == 2

            # Test get_neighbour_coords with sorting
            coords_sorted, verts_sorted = get_neighbour_coords(g, v1, true)
            @test length(coords_sorted) == 2
            @test length(verts_sorted) == 2
        end

        @testset "Weighted degree and geomfilter" begin
            g = DynamicGeometricGraph{2, Float64}()
            coords = [SVector{2, Float64}(10, 10),SVector{2, Float64}(20, 20),
                     SVector{2, Float64}(20, 0),SVector{2, Float64}(0, 25),
                     SVector{2, Float64}(0, 0)]
            for c in coords
                add_vertex!(g, c)
            end
            es = [[1,4],[1,3],[1,2],[1, 5]]

            for v in 1:5
                @test weighted_degree(g, v) == 0
            end

            for (s, d) in es
                add_edge!(g, s, d)
            end

            ds = [weighted_degree(g, v) for v in 1:5]
            @test length(unique(ds)) == 3

            # Test geomfilter
            @test geomfilter(0, 0) < geomfilter(2, 4) < geomfilter(3, 3)
            @test geomfilter(5, 3) == geomfilter(3, 5)  # Symmetric
        end

        @testset "Graph operations" begin
            g = DynamicGeometricGraph{2, Float64}()

            # Test precision
            @test DynamicGeometricGraphs.precision(g) == 1.0

            # Test graph_edits counter
            initial_edits = graph_edits(g)
            add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            @test graph_edits(g) > initial_edits

            # Test freeze (copy)
            v1 = add_vertex!(g, SVector{2, Float64}(1.0, 1.0))
            v2 = add_vertex!(g, SVector{2, Float64}(2.0, 2.0))
            add_edge!(g, v1, v2)

            g2 = freeze(g)
            @test nv(g2) == nv(g)
            @test ne(g2) == ne(g)

            # Modify original, frozen should be unchanged
            v3 = add_vertex!(g, SVector{2, Float64}(3.0, 3.0))
            @test nv(g) != nv(g2)

            # Test Base.copy
            g3 = copy(g)
            @test nv(g3) == nv(g)
        end

        @testset "Clockwise sorting" begin
            # Points around origin
            points = [
                SVector{2, Float64}(1.0, 0.0),   # Right
                SVector{2, Float64}(0.0, 1.0),   # Top
                SVector{2, Float64}(-1.0, 0.0),  # Left
                SVector{2, Float64}(0.0, -1.0)   # Bottom
            ]

            indices = sort_clockwise_indices(points)
            @test length(indices) == 4
            @test all(i in indices for i in 1:4)
        end

        @testset "Custom distance function" begin
            manhattan(a, b) = sum(abs.(a .- b))
            g = DynamicGeometricGraph{2, Float64}(distfun=manhattan, precision=0.5)

            @test DynamicGeometricGraphs.precision(g) == 0.5

            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(1.0, 1.0))
            add_edge!(g, v1, v2)

            # Edge weight should use Manhattan distance
            w = edge_weight(g, v1, v2)
            @test w == 2.0  # |1-0| + |1-0|
        end
    end

    @testset "Bug fixes (GPT-5.1 review 2026-02-17)" begin
        @testset "#4: find_nearest returns correct vertex" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(10.0, 10.0))
            v3 = add_vertex!(g, SVector{2, Float64}(100.0, 100.0))

            # Query point near v1 — must return v1, not last vertex
            query = SVector{2, Float64}(0.1, 0.1)
            @test find_nearest(g, query) == v1

            # Query point near v3
            query2 = SVector{2, Float64}(99.0, 99.0)
            @test find_nearest(g, query2) == v3

            # Exact match should still work
            @test find_nearest(g, SVector{2, Float64}(10.0, 10.0)) == v2

            # Empty graph should error
            g_empty = DynamicGeometricGraph{2, Float64}()
            @test_throws ErrorException find_nearest(g_empty, SVector{2, Float64}(1.0, 1.0))
        end

        @testset "#5: counter uses MVector{2,Int}" begin
            g = DynamicGeometricGraph{2, Float64}()
            @test g.counter isa MVector{2, Int}
            @test g.counter == MVector{2, Int}(0, 0)

            # After operations, counter should still be MVector
            add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            @test g.counter isa MVector{2, Int}

            # Freeze should also produce MVector counter
            g2 = freeze(g)
            @test g2.counter isa MVector{2, Int}
        end

        @testset "#6: translate_graph preserves structure" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(1.0, 0.0))
            v3 = add_vertex!(g, SVector{2, Float64}(0.0, 1.0))
            add_edge!(g, v1, v2)
            add_edge!(g, v2, v3)

            delta = SVector{2, Float64}(5.0, 5.0)
            g2 = translate_graph(g, delta)

            @test nv(g2) == 3
            @test ne(g2) == 2
            # Check coordinates shifted
            crds, vs = get_coords(g2)
            for i in 1:size(crds, 1)
                @test all(crds[i, :] .>= 4.9)  # all coords shifted by 5
            end
        end

        @testset "#6: scale_graph preserves structure" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(1.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(2.0, 0.0))
            add_edge!(g, v1, v2)

            g2 = scale_graph(g, 3.0)
            @test nv(g2) == 2
            @test ne(g2) == 1
            # v1 should be at (3,0), v2 at (6,0)
            crds, vs = get_coords(g2)
            idx1 = findfirst(==(v1), vs)
            @test crds[idx1, 1] ≈ 3.0
        end

        @testset "#6: rotate_graph preserves structure" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(1.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(0.0, 1.0))
            add_edge!(g, v1, v2)

            g2 = rotate_graph(g, Float64(π/2))  # 90 degrees
            @test nv(g2) == 2
            @test ne(g2) == 1
            # (1,0) rotated 90° → (0,1); (0,1) rotated 90° → (-1,0)
            crds, vs = get_coords(g2)
            idx1 = findfirst(==(v1), vs)
            @test crds[idx1, 1] ≈ 0.0 atol=1e-10
            @test crds[idx1, 2] ≈ 1.0 atol=1e-10
        end

        @testset "#7: rem_edge! logic correctness" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(0.0, 0.0))
            v2 = add_vertex!(g, SVector{2, Float64}(1.0, 0.0))
            v3 = add_vertex!(g, SVector{2, Float64}(2.0, 0.0))
            add_edge!(g, v1, v2)

            edits_before = graph_edits(g)

            # Remove existing edge — should increment counter
            rem_edge!(g, v1, v2)
            @test !Graphs.has_edge(g, v1, v2)
            @test graph_edits(g) == edits_before + 1

            edits_after_remove = graph_edits(g)

            # Remove non-existent edge — counter should NOT change
            with_logger(NullLogger()) do
                rem_edge!(g, v1, v2)
            end
            @test graph_edits(g) == edits_after_remove

            # Remove edge with non-existent vertex — should not error
            with_logger(NullLogger()) do
                rem_edge!(g, v1, 999)
            end
            @test graph_edits(g) == edits_after_remove
        end

        @testset "#8: get_vertex_coords/idx return nothing for missing" begin
            g = DynamicGeometricGraph{2, Float64}()
            v1 = add_vertex!(g, SVector{2, Float64}(1.0, 2.0))

            # Valid vertex returns coordinates
            @test get_vertex_coords(g, v1) == SVector{2, Float64}(1.0, 2.0)

            # Invalid vertex returns nothing (not Inf)
            @test get_vertex_coords(g, 999) === nothing

            # Valid coord returns index
            @test get_vertex_idx(g, SVector{2, Float64}(1.0, 2.0)) == v1

            # Invalid coord returns nothing (not Inf)
            @test get_vertex_idx(g, SVector{2, Float64}(99.0, 99.0)) === nothing
        end
    end

end
