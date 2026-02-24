# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Ben Cardoen <b.cardoen@bham.ac.uk>

using Documenter
using DynamicGeometricGraphs

makedocs(
    sitename = "DynamicGeometricGraphs.jl",
    modules = [DynamicGeometricGraphs],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
    authors = "bencardoen",
    repo = "https://github.com/bencardoen/DynamicGeometricGraphs.jl",
    clean = true,
)
