using JSON

nb = JSON.parsefile("Tabellen.ipynb")

# ────────────────────────────────────────────────────────────────────────────
# FIX 1 – Cell 785a100a: round_measurement
#   Änderungen:
#   a) Funktionssignatur: sigdigits::Int=1 als Keyword hinzufügen
#   b) u==0 Guard einfügen (verhindert log10(0) DomainError)
# ────────────────────────────────────────────────────────────────────────────
for cell in nb["cells"]
    if get(cell, "id", "") == "785a100a"
        src = join(cell["source"])
        src = replace(src,
            "using CairoMakie" =>
            "using CairoMakie\nusing LaTeXStrings")
        src = replace(src,
            "function round_measurement(m::Measurement)" =>
            "function round_measurement(m::Measurement; sigdigits::Int=1)")
        src = replace(src,
            "    u = uncertainty(m)\n    if first_significant_digit(u)<3" =>
            "    u = uncertainty(m)\n    # BUGFIX: u==0 würde log10(0) → DomainError\n    if u == 0.0\n        return m\n    end\n    if first_significant_digit(u)<3")
        cell["source"] = [src]
        cell["outputs"] = []
        cell["execution_count"] = nothing
        println("✓ Fixed cell 1 (785a100a): round_measurement")
    end
end

# ────────────────────────────────────────────────────────────────────────────
# FIX 2 – Cell 97f335d7: table()
#   Probleme im Original:
#   a) round_measurement.(sin.(...), sigdigits=3) → MethodError (kein sigdigits)
#   b) sc/bcc/fcc kommen als sortierte Listen zurück, aber danach
#      sc_sorted = [(N, sc[N]) for N in sort(collect(keys(sc)))]
#      behandelt sie als Dict → MethodError
#   c) row1/row2 werden direkt aus Float-Werten erzeugt ohne Measurement-Wrapper
#      → round_measurement bekommt Float statt Measurement → MethodError
# ────────────────────────────────────────────────────────────────────────────
new_table = """function table(th_2theta)
    theta = th_2theta ./ 2
    # Measurement-Wrapper mit Fehler 0
    theta_m = [measurement(Float64(t), 0.0) for t in theta]
    sin_m   = [measurement(sin(Float64(t) * π / 180), 0.0) for t in theta]
    row1 = round_measurement.(theta_m)
    row2 = round_measurement.(sin_m)
    sc, bcc, fcc = generate_reflections(5, length(th_2theta))
    row3 = []
    row4 = Float64[]
    row5 = []
    row6 = Float64[]
    row7 = []
    row8 = Float64[]
    n = min(length(sc), length(bcc), length(fcc))
    for i in 1:n
        push!(row3, sc[i][2][1])
        push!(row4, round(sqrt(Float64(sc[i][1])), sigdigits=3))
        push!(row5, bcc[i][2][1])
        push!(row6, round(sqrt(Float64(bcc[i][1])), sigdigits=3))
        push!(row7, fcc[i][2][1])
        push!(row8, round(sqrt(Float64(fcc[i][1])), sigdigits=3))
    end
    DataFrame(
        theta=value.(row1), sin_theta=value.(row2),
        hkl_sc=row3, sqrtN_sc=row4,
        hkl_bcc=row5, sqrtN_bcc=row6,
        hkl_fcc=row7, sqrtN_fcc=row8
    )
end
"""

for cell in nb["cells"]
    if get(cell, "id", "") == "97f335d7"
        cell["source"] = [new_table]
        cell["outputs"] = []
        cell["execution_count"] = nothing
        println("✓ Fixed cell 5 (97f335d7): table()")
    end
end

# ────────────────────────────────────────────────────────────────────────────
# FIX 3 – Cell 8d7c8d16: table()-Aufruf
#   Probleme im Original:
#   a) 0° in Eingabe → sin(0)=0 → round_measurement(0±0) → log10(0) Crash
#   b) println(th = t[:,2]) → Spaltenzugriff per Index statt Name
#   c) println(root = t[:,8]) → gleicher Fehler + falsche Spalte für Plot
# ────────────────────────────────────────────────────────────────────────────
new_call = """angles_2theta = [30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360]
t = table(angles_2theta)
println(t)
th    = t.sin_theta
root  = t.sqrtN_fcc
delth = 0.0
println("sin(θ) = ", th)
println("√(h²+k²+l²) = ", root)
"""

for cell in nb["cells"]
    if get(cell, "id", "") == "8d7c8d16"
        cell["source"] = [new_call]
        cell["outputs"] = []
        cell["execution_count"] = nothing
        println("✓ Fixed cell 6 (8d7c8d16): table() Aufruf")
    end
end

# ────────────────────────────────────────────────────────────────────────────
# FIX 4 – Cell 94348c22: letzter Plot – Skalierung reparieren
#   Probleme im Original:
#   a) nice_ticks_and_labels schlägt fehl wenn alle y-Werte = 0 (N=0 Reflex)
#   b) alignedlabels gibt leere Strings → Achse unbeschriftet
#   c) x=th, y=root: root enthält 0.0 (N=0) → log10(0) crash in nice_ticks
#   Lösung: Ersten Eintrag (N=0, hkl=(0,0,0)) überspringen,
#            Skalierung direkt über xlims!/ylims! setzen
# ────────────────────────────────────────────────────────────────────────────
new_plot = raw"""function nice_ticks_and_labels(x, n_labels=11, sigdigit=8)
    xmin, xmax = minimum(x), maximum(x)
    if xmin == xmax
        return [xmin], [string(xmin)]
    end
    order = Float64(floor(Int, log10(xmax - xmin + 1e-12)))
    n = Float64(floor(Int, xmax / 10^order)) + 1
    tick_max = n * 10^order
    ticks  = round.(collect(range(0.0, stop=tick_max, step=5*10^(order-2))); digits=sigdigit)
    labels = collect(range(0.0, stop=tick_max, length=n_labels))
    return ticks, labels
end

function alignedlabels(ticks, labels; atol=0.06)
    aligned = String[]
    for x in ticks
        match = findfirst(lbl -> isapprox(x, lbl; atol=atol, rtol=0), labels)
        if match === nothing
            push!(aligned, "")
        else
            push!(aligned, string(labels[match]))
        end
    end
    return aligned
end

# Ersten Eintrag (N=0, hkl=(0,0,0)) weglassen – physikalisch kein Bragg-Reflex
x = th[2:end]       # sin(θ)
y = root[2:end]     # √(h²+k²+l²)
Δx = delth          # Fehler auf sin(θ)

with_theme(theme_latexfonts()) do
    fig = Figure(size=(700, 500))
    ax  = Axis(fig[1,1],
        title  = L"Bragg-Auswertung: $\\sqrt{N}$ vs $\\sin\\theta$",
        xlabel = L"$\\sin(\\theta)$",
        ylabel = L"$\\sqrt{h^2+k^2+l^2}$")

    # Skalierung automatisch aus Datenpunkten
    xpad = 0.05 * (maximum(x) - minimum(x) + 1e-6)
    ypad = 0.1  * (maximum(y) - minimum(y) + 1e-6)
    xlims!(ax, minimum(x) - xpad, maximum(x) + xpad)
    ylims!(ax, 0.0, maximum(y) + ypad)

    # Fehlerbalken nur in x (Δx = delth)
    errorbars!(ax, x, y, fill(Δx, length(x)), color=:grey, whiskerwidth=5, direction=:x)

    # Datenpunkte
    scatter!(ax, x, y, markersize=8, color=:black)

    fig
end
"""

for cell in nb["cells"]
    if get(cell, "id", "") == "94348c22"
        cell["source"] = [new_plot]
        cell["outputs"] = []
        cell["execution_count"] = nothing
        println("✓ Fixed cell 7 (94348c22): Plot-Skalierung")
    end
end

# Notebook speichern
open("Tabellen.ipynb", "w") do f
    JSON.print(f, nb, 1)
end
println("\n✅ Tabellen.ipynb erfolgreich gespeichert!")
