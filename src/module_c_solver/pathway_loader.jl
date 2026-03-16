using JSON
using DataFrames
using JuMP

"""
Module C: Solver Handshake
Optimized for Julia 1.12+ with robust pathing.
"""
function load_pathway_data(filepath::String)
    # Check if file exists before opening to provide better error messages
    if !isfile(filepath)
        error("FileNotFound: Could not find JSON at $filepath. Ensure Module B ran successfully.")
    end
    
    raw_data = JSON.parsefile(filepath)
    nodes = raw_data["nodes"]
    edges = raw_data["edges"]

    println("[*] Ingested $(length(nodes)) nodes.")

    df_nodes = DataFrame(
        id = [n["id"] for n in nodes],
        name = [n["name"] for n in nodes],
        type = [n["type"] for n in nodes]
    )

    return df_nodes, edges
end

# --- Robust Execution Block ---
function main()
    # @__DIR__ points to src/module_c_solver/
    # We go up one level to src/, then up to root, then into data/
    json_path = joinpath(@__DIR__, "..", "..", "data", "enriched_pathway.json")
    
    pathway_df, edges = load_pathway_data(json_path)

    model = Model()
    @variable(model, 0.01 <= activity[pathway_df.id] <= 100.0)

    println("[✔] Julia Handshake Complete. $(length(activity)) variables ready.")
end

# Wrap the execution in invokelatest to satisfy Julia 1.12 world age semantics
Base.invokelatest(main)