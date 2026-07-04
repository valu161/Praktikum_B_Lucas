using JSON
using Test

# Load notebook
nb = JSON.parsefile("Tabellen.ipynb")

println("Executing notebook cells...")
for (i, cell) in enumerate(nb["cells"])
    if cell["cell_type"] == "code"
        id = get(cell, "id", "unknown")
        src = join(cell["source"])
        if isempty(strip(src))
            continue
        end
        println("\n--- Executing Cell $i (id: $id) ---")
        try
            # Evaluate code in Main module using include_string
            include_string(Main, src, "Cell-$i")
            println("✓ Cell $i executed successfully")
        catch e
            println("❌ Error in Cell $i:")
            showerror(stdout, e)
            println()
            # Don't rethrow, let's see how far we get or exit
            rethrow(e)
        end
    end
end
println("\nAll cells executed successfully!")
