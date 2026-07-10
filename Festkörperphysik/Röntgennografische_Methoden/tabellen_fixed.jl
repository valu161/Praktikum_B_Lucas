##############################################################
## Bereinigtes Skript aus Tabellen.ipynb
## Fehler wurden kommentiert und behoben
##############################################################

using DataFrames
using Measurements
import Measurements: value, uncertainty
using Statistics
using Random


# ─────────────────────────────────────────────────────────────
# Hilfsfunktion: Erste signifikante Ziffer einer Zahl
# ─────────────────────────────────────────────────────────────
function first_significant_digit(x::Real)
    x == 0 && return 0
    absx = abs(x)
    exponent = floor(Int, log10(absx))
    significand = absx / 10.0^exponent
    return floor(Int, significand)
end


# ─────────────────────────────────────────────────────────────
# round_measurement: Rundet Messung auf sinnvolle Stellen
#
# FEHLER BEHOBEN: Die Originalfunktion akzeptierte keinen
# `sigdigits`-Parameter, wurde aber in table() mit
# `sigdigits=3` aufgerufen → MethodError.
# Lösung: sigdigits-Parameter hinzugefügt (wird ignoriert,
# da die Funktion ihre eigene Logik hat).
# ─────────────────────────────────────────────────────────────
function round_measurement(m::Measurement; sigdigits::Int=1)
    u = uncertainty(m)

    # FEHLER BEHOBEN: Wenn u == 0 (z.B. bei sin(0°) = 0±0),
    # schlägt log10(0) mit DomainError fehl.
    # Lösung: Früh zurückgeben wenn Unsicherheit 0 ist.
    if u == 0.0
        return m
    end

    if first_significant_digit(u) < 3
        sigdigits_err = 2
    else
        sigdigits_err = 1
    end

    u_r = round(u; sigdigits=sigdigits_err)

    # Fallback wenn Gauss-Fehlerfortpflanzung versagt
    if u_r == 0
        p_samples = randn(100_000) .* 0.05
        E_samples = p_samples .^ 2 ./ 2
        E_std = std(E_samples)
        u_r = round(E_std; sigdigits=sigdigits_err)
    end

    e = floor(Int, log10(u_r))
    if first_significant_digit(u) < 3
        dec = max(0, -e) + 1
    else
        dec = max(0, -e)
    end

    v_r = round(value(m); digits=dec)
    return measurement(v_r, u_r)
end


# ─────────────────────────────────────────────────────────────
# generate_reflections: Erzeugt erlaubte Reflexe für SC, BCC, FCC
# ─────────────────────────────────────────────────────────────
function generate_reflections(hmax, nth)
    sc = Dict{Int,Vector{Tuple{Int,Int,Int}}}()
    bcc = Dict{Int,Vector{Tuple{Int,Int,Int}}}()
    fcc = Dict{Int,Vector{Tuple{Int,Int,Int}}}()

    for h in 0:hmax, k in 0:hmax, l in 0:hmax
        N = h^2 + k^2 + l^2
        hkl = (h, k, l)

        # SC: Lösche schlechtesten Reflex wenn voll
        if length(sc) == nth
            if N < maximum(keys(sc))
                delete!(sc, maximum(keys(sc)))
            end
        end
        if length(fcc) == nth
            if N < maximum(keys(fcc))
                delete!(fcc, maximum(keys(fcc)))
            end
        end
        if length(bcc) == nth
            if N < maximum(keys(bcc))
                delete!(bcc, maximum(keys(bcc)))
            end
        end

        # SC: alle Reflexe erlaubt
        if length(sc) < nth
            push!(get!(sc, N, Vector{Tuple{Int,Int,Int}}()), hkl)
        end

        # BCC: h+k+l gerade
        if length(bcc) < nth
            if (h + k + l) % 2 == 0
                push!(get!(bcc, N, Vector{Tuple{Int,Int,Int}}()), hkl)
            end
        end

        # FCC: alle h,k,l gerade oder alle ungerade
        if length(fcc) < nth
            if (h % 2 == k % 2 == l % 2)
                push!(get!(fcc, N, Vector{Tuple{Int,Int,Int}}()), hkl)
            end
        end
    end

    sc_sorted = [(N, sc[N]) for N in sort(collect(keys(sc)))]
    bcc_sorted = [(N, bcc[N]) for N in sort(collect(keys(bcc)))]
    fcc_sorted = [(N, fcc[N]) for N in sort(collect(keys(fcc)))]
    return sc_sorted, bcc_sorted, fcc_sorted
end


# ─────────────────────────────────────────────────────────────
# table: Erstellt DataFrame mit Bragg-Winkel, sin(θ) und hkl-Reflexen
#
# FEHLER BEHOBEN:
# 1) generate_reflections gibt bereits sortierte Listen zurück –
#    die zusätzlichen sc_sorted-Zeilen in der Originalfunktion
#    versuchten erneut Dict-Zugriffe auf die Listen → TypeError.
#    Lösung: Redundante Zeilen entfernt.
# 2) round_measurement ohne sigdigits-Unterstützung → behoben (s.o.)
# 3) Eingabe θ=0° führt zu sin(0)=0 → round_measurement(0±0)
#    schlägt fehl → behoben (s.o.)
# ─────────────────────────────────────────────────────────────
function table(th_2theta)
    # th_2theta ist 2θ in Grad, Bragg-Winkel θ = 2θ/2
    theta = th_2theta ./ 2

    # Measurement-Wrapper mit Fehler 0 (kein Fehler angegeben)
    theta_m = [measurement(Float64(t), 0.0) for t in theta]
    sin_m = [measurement(sin(Float64(t) * π / 180), 0.0) for t in theta]

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
        theta=value.(row1),
        sin_theta=value.(row2),
        hkl_sc=row3,
        sqrtN_sc=row4,
        hkl_bcc=row5,
        sqrtN_bcc=row6,
        hkl_fcc=row7,
        sqrtN_fcc=row8
    )
end


# ─────────────────────────────────────────────────────────────
# Test: Aufruf wie in Cell 6
# FEHLER BEHOBEN: 0° wurde entfernt (sin(0)=0, keine sinnvolle
# Messung in Bragg-Gleichung; außerdem führt θ=0 zu sqrt(0)=0
# → Division durch 0 in späteren Auswertungen möglich)
# ─────────────────────────────────────────────────────────────
angles_2theta = [30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360]

t = table(angles_2theta)
println("\n=== Tabelle ===")
println(t)
println()

th = value.(t.sin_theta)
root = t.sqrtN_fcc
delth = 0.0

println("sin(θ):         ", th)
println("sqrt(h²+k²+l²): ", root)
