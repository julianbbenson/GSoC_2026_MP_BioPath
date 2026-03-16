using JSON
using GraphRecipes
using Plots

# Use the stable GR backend
gr()

function render_pathway(json_path::String)
    data = JSON.parsefile(json_path)
    nodes = data["nodes"]
    edges = data["edges"]

    # Map IDs to Indices
    node_ids = [n["id"] for n in nodes]
    id_to_idx = Dict(id => i for (i, id) in enumerate(node_ids))
    n = length(nodes)
    adj = zeros(Int, n, n)
    for e in edges
        s, d = get(id_to_idx, e["source"], 0), get(id_to_idx, e["target"], 0)
        if s > 0 && d > 0 adj[s, d] = 1 end
    end

    println("[*] Rendering final stable constellation...")

    # AESTHETIC FIXES:
    # 1. Removed edgestyle=:curved (unsupported by GR)
    # 2. Kept markercolor and markersize for legibility
    # 3. Size is optimized for 69 nodes without crashing the buffer
    p = graphplot(adj, 
                  names=[n["name"] for n in nodes], 
                  method=:spring,
                  markercolor=:lightcyan,
                  markerstrokecolor=:deepskyblue,
                  markerstrokewidth=1.2,
                  markershape=:rect,
                  markersize=0.12,
                  linecolor=:gray,
                  linealpha=0.3,
                  fontsize=6, 
                  size=(1600, 1600), 
                  title="MP-BioPath Stitched Network (R-HSA-5673001)")

    out_path = joinpath(@__DIR__, "..", "..", "docs", "stitched_constellation_final.png")
    
    # Save directly; the crash was the cause of the empty file, not the timing.
    savefig(p, out_path)
    
    println("[✔] Final constellation saved to: $out_path")
end

function main()
    json_path = joinpath(@__DIR__, "..", "..", "data", "enriched_pathway.json")
    render_pathway(json_path)
end

Base.invokelatest(main)