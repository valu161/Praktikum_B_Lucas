#Diesen Code vor Gaußfits ausführen. Mit linksklick die peks markieren und mit 
#maus fenster schließen. Der Code ist von Julian, da ich zu kacke bin das Zeug selber 
#ordentlich zu mache






using GLMakie
using DelimitedFiles
using Printf # For formatting file names
using JSON # Für das Speichern der Peak-Positionen in einer Datei
using Statistics # Für mean im späteren Fit-Teil
GLMakie.activate!()
# --- 1. Function to read .xye data ---
function read_xye_data(filepath::String)
    if !isfile(filepath)
        @error "File not found at $filepath"
        return nothing, nothing
    end

    lines = readlines(filepath)
    data_lines = lines[2:end]

    two_theta_values = Float64[]
    intensity_values = Float64[]

    for line in data_lines
        parts = split(strip(line))
        if length(parts) == 3
            try
                push!(two_theta_values, parse(Float64, parts[1]))
                push!(intensity_values, parse(Float64, parts[2]))
            catch e
                @warn "Could not parse line: $line ($e)"
                continue
            end
        end
    end
    return two_theta_values, intensity_values
end

# --- 2. File paths for your measurements ---
file_paths = Dict(
    "Probe 1" => "Messdaten_Pulverdiffraktonomie/Probe1_Messdaten.xye",
    "Probe 2" => "Messdaten_Pulverdiffraktonomie/Probe2_Messdaten_neu.xye",
    "Probe 3" => "Messdaten_Pulverdiffraktonomie/Probe3_Messdaten.xye",
)

# --- 3. Interactive Plotting Function ---
function interactive_peak_selection(probe_name::String, filepath::String, all_probes_peak_data::Dict)
    x, y = read_xye_data(filepath)

    if x === nothing || y === nothing || isempty(x)
        println("Skipping interactive selection for $probe_name due to missing or invalid data.")
        return
    end

    f = Figure(size=(1000, 600))
    ax = Axis(f[1, 1],
        xlabel="2θ (°)",
        ylabel="Intensität (log)",
        title="Klicke auf Peaks für $probe_name (Linksklick, ESC zum Schließen)",
        yscale=log10 # Logarithmische Skala auf der y-Achse
    )

    lines!(ax, x, y, linewidth=1.5, color=:blue)

    # --- Korrekturen hier ---
    current_probe_peaks = Float64[]
    plotted_markers = Observable(Point2f[])

    scatter!(ax, plotted_markers, marker='⬇', markersize=20, color=:red, strokewidth=1, strokecolor=:black)

    # Binde den Event-Handler direkt an die Achse (ax.scene) statt an die Figur (f)
    on(events(ax.scene).mousebutton) do event
        # Nur auf den linken Mausklick reagieren
        if event.button == Mouse.left && event.action == Mouse.press
            # mouseposition(ax.scene) gibt die Position direkt in Datenkoordinaten zurück!
            # Das macht die Konvertierung von Pixel zu Daten überflüssig.
            clicked_2theta, _ = mouseposition(ax.scene)

            # Füge den Peak zur Liste hinzu
            push!(current_probe_peaks, clicked_2theta)
            @info @sprintf("Geklickte 2θ-Position für %s: %.4f°", probe_name, clicked_2theta)

            # --- Visuelles Feedback ---
            # Vertikale Linie an der geklickten Position
            vlines!(ax, clicked_2theta, color=:red, linestyle=:dot, linewidth=1.5)

            # Marker an der Oberseite des Plots hinzufügen (modernisierter Ansatz)
            # 1. Hole die aktuelle obere Y-Grenze der Achse
            # === GEÄNDERTE STELLE START ===
            ylims = ax.finallimits[] # Nutze finallimits[] statt limits.val, da limits.val (nothing, nothing) sein kann und keine origin/widths besitzt
            current_y_limit_top = ylims.origin[2] + ylims.widths[2]
            # === GEÄNDERTE STELLE ENDE ===

            # 2. Erstelle den neuen Punkt
            new_marker = Point2f(clicked_2theta, current_y_limit_top * 0.95)

            # 3. Aktualisiere das Observable auf die empfohlene Weise
            # Dies fügt den neuen Marker zum bestehenden Array hinzu und benachrichtigt den Plot
            # === GEÄNDERTE STELLE START ===
            push!(plotted_markers[], new_marker)
            notify(plotted_markers) # Nutze notify() um das Observable sauber zu aktualisieren
            # === GEÄNDERTE STELLE ENDE ===
        end
        # Gib `Consume(false)` zurück, damit andere Interaktionen (z.B. Zoomen) weiterhin funktionieren
        return Consume(false)
    end

    println("Interaktives Fenster für $probe_name geöffnet. Klicke auf die Peak-Positionen.")
    println("Drücke ESC, um das Fenster zu schließen und fortzufahren.")

    display(f)
    wait(f.scene) # Warten, bis das Fenster geschlossen wird
    println("Fenster für $probe_name geschlossen. Gesammelte Peaks: $(current_probe_peaks)")

    # Speichern der gesammelten Peaks für diese Probe
    all_probes_peak_data[probe_name] = current_probe_peaks
end

# --- Main execution block ---
function main_interactive_selection()
    global_peak_data = Dict{String,Vector{Float64}}() # Dictionary zum Speichern aller Peak-Daten

    # === GEÄNDERTE STELLE START ===
    # Sortiere die Schlüssel alphabetisch, damit die Proben der Reihe nach (Probe 1, Probe 2, Probe 3) geöffnet werden
    for probe_name in sort(collect(keys(file_paths)))
        filepath = file_paths[probe_name]
        interactive_peak_selection(probe_name, filepath, global_peak_data)
    end
    # === GEÄNDERTE STELLE ENDE ===

    # Speichere die gesammelten Peak-Positionen in einer JSON-Datei
    output_json_path = "peak_positions.json"
    open(output_json_path, "w") do f
        JSON.print(f, global_peak_data, 2) # Pretty print with indent 2
    end
    println("\nAlle Peak-Positionen wurden in '$output_json_path' gespeichert.")
    println("Nutze diese Datei, um deine `peak_configs` für den Fit-Code zu aktualisieren.")
    println("Interaktive Peak-Auswahl für alle Proben abgeschlossen.")
end

# Führe die interaktive Auswahl aus
main_interactive_selection()