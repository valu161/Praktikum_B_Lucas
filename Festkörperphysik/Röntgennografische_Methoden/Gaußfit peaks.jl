#Erst in Get peaks den Code ausführen uns im Interactiven Fenster die peaks 
#markieren. Dann den code hier ausführen.


using LsqFit
using GLMakie
using CairoMakie      # Für statische Plots und das Speichern als PDF
using DelimitedFiles
using Printf
using JSON
using Statistics
using CSV             # Zum Speichern von Daten in CSV-Dateien
using DataFrames      # Zum einfachen Arbeiten mit tabellarischen Daten für CSV
CairoMakie.activate!()      # jetzt ist CairoMakie das aktive Backend
set_theme!(fontsize=25)

# --- 1. Hilfsfunktion zum Einlesen der .xye Daten ---
function read_xye_data(filepath::String)
    if !isfile(filepath)
        @error "Datei nicht gefunden unter $filepath"
        return nothing, nothing
    end

    lines = readlines(filepath)
    data_lines = lines[2:end] # Erste Zeile ist der Header

    two_theta_values = Float64[]
    intensity_values = Float64[]

    for line in data_lines
        parts = split(strip(line))
        if length(parts) >= 2
            try
                push!(two_theta_values, parse(Float64, parts[1]))
                push!(intensity_values, parse(Float64, parts[2]))
            catch e
                @warn "Konnte Zeile nicht parsen: $line ($e)"
                continue
            end
        else
            @warn "Ungültiges Zeilenformat übersprungen (zu wenige Spalten): $line"
        end
    end
    return two_theta_values, intensity_values
end

# --- 2. Definition der Modellfunktionen ---
# Eine einzelne Gaußkurve
gaussian(x, A, mu, sigma) = A * exp.(-(x .- mu) .^ 2 / (2 * sigma^2))

# Lineare Basislinie
linear_baseline(x, a, b) = a .* x .+ b

# Modellfunktion für eine Gruppe von Gaußkurven MIT BASISLINIE (für den lokalen Fit)
# p: Parameter für diese spezifische Gruppe von Gaußkurven UND Basislinie
# num_gaussians: Wie viele Gaußkurven in dieser Gruppe sind
function local_peak_model_with_baseline(x_data, p, num_gaussians::Int)
    # BEHOBEN: Nutzt den Typ von p (z.B. ForwardDiff.Dual), statt fest Float64 vorzugeben
    y_model = zeros(eltype(p), length(x_data))

    # Die letzten zwei Parameter sind für die Basislinie: a (Steigung), b (Achsenabschnitt)
    idx_baseline_a = length(p) - 1
    idx_baseline_b = length(p)

    a_b = p[idx_baseline_a]
    b_b = p[idx_baseline_b]

    y_model .+= linear_baseline(x_data, a_b, b_b)

    param_idx = 1 # Startet bei den Gauß-Parametern
    for i in 1:num_gaussians
        # Sicherstellen, dass genügend Parameter für die Gaußkurve vorhanden sind
        if param_idx + 2 > length(p) - 2
            @warn "Nicht genügend Gauß-Parameter in 'p' für erwartete Gaußkurve $i."
            break
        end
        A, mu, sigma = p[param_idx], p[param_idx+1], p[param_idx+2]
        y_model .+= gaussian(x_data, A, mu, sigma)
        param_idx += 3
    end
    return y_model
end

# NEUE FUNKTION: Speichert nur den gewichteten Mittelwert und dessen Fehler
function save_weighted_mean_results_to_csv(probe_name::String, results::Vector{<:NamedTuple})
    output_dir = "Fittergebnisse_CSV"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    filename = @sprintf("%s/%s_Gewichteter_Mittelwert.csv", output_dir, replace(probe_name, " " => "_"))

    df = DataFrame(
        Peak_Label=[r.peak_label for r in results],
        Weighted_Mean_Theta=[r.weighted_mean_theta for r in results],
        Error_Weighted_Mean_Theta=[r.error_weighted_mean_theta for r in results]
    )

    try
        CSV.write(filename, df)
        println("Gewichtete Mittelwerte für $probe_name in $filename gespeichert.")
    catch e
        @error "Fehler beim Speichern der CSV-Datei für $probe_name: $e"
    end
end


# --- 3. Hauptfunktion für Fitting und Plotting (KORRIGIERT) ---
function fit_and_plot_peaks(probe_name::String, filepath::String, json_peak_positions::AbstractDict)
    x_obs, y_obs = read_xye_data(filepath)

    if x_obs === nothing || y_obs === nothing || isempty(x_obs)
        println("Überspringe Fitting für $probe_name: Keine gültigen Daten gefunden.")
        return
    end

    # Lese die Peak-Positionen aus dem JSON (hier sind es wieder nur die mus)
    initial_mu_guesses = get(json_peak_positions, probe_name, Float64[])
    if isempty(initial_mu_guesses)
        @warn "Keine Peak-Positionen für '$probe_name' in 'peak_positions.json' gefunden. Kann keinen Fit durchführen."
        return
    end

    peak_configs = []
    sort!(initial_mu_guesses)

    # *** NEU: Hier kannst du den Breitenfaktor anpassen! ***
    PEAK_RANGE_WIDTH_FACTOR = 1 # Passe diesen Wert an, z.B. 1.0, 2.0 etc.
    DEFAULT_NUM_GAUSSIANS = 2 # Standardanzahl der Gaußkurven pro Peak

    for mu_guess in initial_mu_guesses
        search_range_min = mu_guess - 0.7
        search_range_max = mu_guess + 0.7

        idx_in_range = findall(x -> search_range_min < x < search_range_max, x_obs)

        initial_amplitude = if isempty(idx_in_range)
            maximum(y_obs) / 10
        else
            y_obs[idx_in_range[argmax(y_obs[idx_in_range])]] * 1.05
        end
        initial_amplitude = max(initial_amplitude, 1.0)

        initial_baseline_b = if isempty(idx_in_range)
            clamp(mean(y_obs) * 0.01, 1e-3, Inf)
        else
            clamp(mean(y_obs[idx_in_range]) * 0.1, 1e-3, Inf)
        end

        initial_baseline_a = 0.0

        current_peak_range = [mu_guess - PEAK_RANGE_WIDTH_FACTOR, mu_guess + PEAK_RANGE_WIDTH_FACTOR]

        num_gaussians_for_current_peak = DEFAULT_NUM_GAUSSIANS
        if abs(mu_guess - 38.02) < 0.1 # Beispiel: Spezifisch für diesen Peak
            num_gaussians_for_current_peak = 1
        end

        # Erstelle die initialen Parameter für diesen Peak
        initial_params = Vector{Float64}(undef, 3 * num_gaussians_for_current_peak + 2)

        for i in 1:num_gaussians_for_current_peak
            initial_params[3*(i-1)+1] = initial_amplitude / (i == 1 ? 1.0 : 1.5)
            initial_params[3*(i-1)+2] = mu_guess + (i - 1) * 0.05
            initial_params[3*(i-1)+3] = 0.1
        end

        initial_params[end-1] = initial_baseline_a
        initial_params[end] = initial_baseline_b

        push!(peak_configs, (
            range_2theta=current_peak_range,
            num_gaussians=num_gaussians_for_current_peak,
            initial_params=initial_params,
            label=@sprintf("Peak bei %.2f°", mu_guess)
        ))
    end

    fitted_baselines_and_params = Dict{String,Any}()
    weighted_mean_theta_results = Vector{NamedTuple{(:peak_label, :weighted_mean_theta, :error_weighted_mean_theta),Tuple{String,Float64,Float64}}}()


    println("\n### Bearbeite Probe: $probe_name ###")
    println("--- Führe lokale Fits durch ---")

    for (peak_idx, config) in enumerate(peak_configs)
        range_start, range_end = config.range_2theta
        idx_local_range = findall(x -> range_start <= x <= range_end, x_obs)
        x_local = x_obs[idx_local_range]
        y_local = y_obs[idx_local_range]

        if isempty(x_local)
            @warn "       Keine Datenpunkte im Bereich $(range_start)° - $(range_end)° für $(config.label). Überspringe diesen Fit."
            continue
        end

        local_model_func(x_data, p) = local_peak_model_with_baseline(x_data, p, config.num_gaussians)

        try
            local_fit = curve_fit(local_model_func, x_local, y_local, config.initial_params, maxIter=500)
            local_fitted_params = local_fit.param

            # --- Standardfehler aus der Kovarianzmatrix extrahieren ---
            param_errors = try
                stderror(local_fit)
            catch e
                @warn "       Konnte Standardfehler für $(config.label) nicht berechnen: $e. Setze Fehler auf NaN."
                fill(NaN, length(local_fitted_params))
            end

            param_idx_current_local = 1

            # Listen für unterschiedliche Zwecke:
            current_gauss_params_for_weighted_mean = [] # Für die Berechnung des gewichteten Mittelwerts und dessen Fehler
            current_gauss_params_for_plotting = []      # Für die Rekonstruktion der Gesamt-Fit-Linie im Plot

            for i in 1:config.num_gaussians
                if param_idx_current_local + 2 > length(local_fitted_params) - 2
                    @warn "       Nicht genügend gefittete Gauß-Parameter für Komponente $i von $(config.label). Möglicher Fit-Fehler."
                    break
                end

                A_fit = local_fitted_params[param_idx_current_local]
                mu_fit = local_fitted_params[param_idx_current_local+1]
                sigma_fit = local_fitted_params[param_idx_current_local+2]

                err_A = param_errors[param_idx_current_local]
                err_mu = param_errors[param_idx_current_local+1]

                if A_fit > 0 && sigma_fit > 1e-6 && !isnan(mu_fit) && !isnan(err_mu) && !isinf(err_mu)
                    # Sammeln für die gewichtete Mittelwertsberechnung
                    push!(current_gauss_params_for_weighted_mean, (mu=mu_fit, A=A_fit, err_mu=err_mu))
                    # Sammeln für das Plotten der Fit-Linie
                    push!(current_gauss_params_for_plotting, (A=A_fit, mu=mu_fit, sigma=sigma_fit))
                else
                    @warn "       Gauß-Komponente $i von $(config.label) hatte unrealistische Fit-Parameter oder Fehler. Wird übersprungen."
                end
                param_idx_current_local += 3
            end

            # Speichern der Fit-Parameter für den Plot
            if !isempty(current_gauss_params_for_plotting)
                a_baseline_fit = local_fitted_params[end-1]
                b_baseline_fit = local_fitted_params[end]
                fitted_baselines_and_params[config.label] = (a_baseline_fit, b_baseline_fit, current_gauss_params_for_plotting, range_start, range_end)
            end

            # Berechnung des gewichteten Mittelwerts und dessen Fehler für diesen Peak
            if !isempty(current_gauss_params_for_weighted_mean)
                if length(current_gauss_params_for_weighted_mean) == 1
                    weighted_mean_2theta = current_gauss_params_for_weighted_mean[1].mu
                    error_weighted_mean_2theta = current_gauss_params_for_weighted_mean[1].err_mu
                else
                    weights = [g.A for g in current_gauss_params_for_weighted_mean]
                    mus_2theta = [g.mu for g in current_gauss_params_for_weighted_mean]
                    err_mus_2theta = [g.err_mu for g in current_gauss_params_for_weighted_mean]

                    if sum(weights) < 1e-9
                        @warn "       Summe der Amplituden für $(config.label) ist zu klein. Gewichteter Mittelwert wird nicht berechnet."
                        continue
                    end

                    weighted_mean_2theta = sum(mus_2theta .* weights) / sum(weights)

                    sum_sq_weighted_errors = sum((w^2 * err_mu^2) for (w, err_mu) in zip(weights, err_mus_2theta))
                    error_weighted_mean_2theta = sqrt(sum_sq_weighted_errors) / sum(weights)
                end

                # Umrechnung nach Theta und Fehleranpassung
                weighted_mean_theta = weighted_mean_2theta / 2.0
                error_weighted_mean_theta = error_weighted_mean_2theta / 2.0

                push!(weighted_mean_theta_results, (
                    peak_label=config.label,
                    weighted_mean_theta=weighted_mean_theta,
                    error_weighted_mean_theta=error_weighted_mean_theta
                ))
                println("       $(config.label): Gewichteter Theta = $(@sprintf("%.4f", weighted_mean_theta))° ± $(@sprintf("%.4f", error_weighted_mean_theta))°")
            end

        catch e
            @error "    Fit für $(config.label) FEHLGESCHLAGEN: $e"
            println("    Bitte initial_params oder range_2theta für diesen Peak überprüfen. Dieser Peak wird nicht im Plot berücksichtigt.")
        end
    end

    # --- Plotting ---
    f = Figure(size=(1000, 700))

    min_y_data = minimum(y_obs)
    max_y_data = maximum(y_obs)

    plot_ymin = max(0.1, min_y_data * 0.9)
    plot_ymax = max_y_data * 2

    ax = Axis(f[1, 1],
        xlabel="2θ (°)",
        ylabel="Intensität (log)",
        title="Diffraktogramm und Fit: $probe_name",
        yscale=log10,
        limits=(minimum(x_obs), maximum(x_obs), plot_ymin, plot_ymax)
    )

    y_obs_clipped = clamp.(y_obs, 0.1, Inf)
    lines!(ax, x_obs, y_obs_clipped, linewidth=1.0, color=:grey, label="Original Daten")

    total_fitted_y_model = fill(NaN, length(x_obs))

    # Plotten der roten Fit-Linie
    for (label, fit_data) in fitted_baselines_and_params
        a_b, b_b, gauss_params_for_plot, range_start, range_end = fit_data # Hier verwenden wir die Parameter für den Plot
        idx_local_range = findall(x -> range_start <= x <= range_end, x_obs)
        x_local = x_obs[idx_local_range]

        if !isempty(x_local)
            local_model_values = linear_baseline(x_local, a_b, b_b)
            for g in gauss_params_for_plot # Jede Gauß-Komponente für diesen Peak
                local_model_values .+= gaussian(x_local, g.A, g.mu, g.sigma)
            end
            total_fitted_y_model[idx_local_range] = local_model_values
        end
    end

    total_fitted_y_model_clipped = clamp.(total_fitted_y_model, 0.1, Inf)
    lines!(ax, x_obs, total_fitted_y_model_clipped, linestyle=:dash, linewidth=2.0, color=:red, label="Gesamt-Fit")

    # Plotten der berechneten gewichteten Mittelwerte im Diagramm
    if !isempty(weighted_mean_theta_results)
        final_limits = ax.finallimits[]
        plot_ymin_actual = final_limits.origin[2]
        plot_ymax_actual = final_limits.origin[2] + final_limits.widths[2]

        for res in weighted_mean_theta_results
            weighted_mean_2theta = res.weighted_mean_theta * 2          # θ  → 2θ

            # Index des Datenpunkts, der der 2θ‑Position am nächsten liegt
            _, idx_nearest = findmin(abs.(x_obs .- weighted_mean_2theta))

            # Höhe des Fits (alternativ: y_obs_clipped[idx_nearest])
            y_peak = total_fitted_y_model_clipped[idx_nearest]
            text_y_pos = clamp(y_peak * 1.15,                           # 15 % darüber
                plot_ymin_actual * 1.05,
                plot_ymax_actual * 0.95)

            text!(ax,
                @sprintf("%.2f°", weighted_mean_2theta),
                position=(weighted_mean_2theta, text_y_pos),
                align=(:center, :bottom),
                color=:black,
                fontsize=12)
        end
    end


    axislegend(ax, [lines!(ax, [0], [0], color=:grey), lines!(ax, [0], [0], color=:red)], ["Original Daten", "Gesamt-Fit"], position=:rt)

    output_dir_pdfs = "Gefittete_Diffraktogramme_PDFs"
    if !isdir(output_dir_pdfs)
        mkpath(output_dir_pdfs)
    end
    filename_pdf = @sprintf("%s/%s_Fitted_Diffraktogramm.pdf", output_dir_pdfs, replace(probe_name, " " => "_"))

    try
        CairoMakie.save(filename_pdf, f)
        println("Plot für $probe_name in $filename_pdf gespeichert.")
    catch e
        @error "Fehler beim Speichern des Plots für $probe_name: $e"
        println("Bitte stellen Sie sicher, dass CairoMakie korrekt installiert und geladen ist.")
    end

    # --- Aufruf der spezifischen Speicherfunktion für den gewichteten Mittelwert ---
    if !isempty(weighted_mean_theta_results)
        save_weighted_mean_results_to_csv(probe_name, weighted_mean_theta_results)
    else
        println("Keine gewichteten Mittelwerte zum Speichern in CSV für $probe_name vorhanden.")
    end
end

# --- Hauptausführung für alle Proben ---
function main_fitting_process()
    json_peak_path = "peak_positions.json"
    if !isfile(json_peak_path)
        @error "Fehler: 'peak_positions.json' wurde nicht gefunden. Bitte stelle sicher, dass die Datei existiert und korrekt formatiert ist."
        return
    end

    json_peak_data = try
        open(json_peak_path, "r") do f
            JSON.parse(f)
        end
    catch e
        @error "Fehler beim Laden oder Parsen von 'peak_positions.json': $e. Bitte stelle sicher, dass die Datei im JSON-Format ist."
        return
    end

    file_paths = Dict(
        "Probe 1" => "Messdaten_Pulverdiffraktonomie/Probe1_Messdaten.xye",
        "Probe 2" => "Messdaten_Pulverdiffraktonomie/Probe2_Messdaten_neu.xye",
        "Probe 3" => "Messdaten_Pulverdiffraktonomie/Probe3_Messdaten.xye",
    )

    # === GEÄNDERTE STELLE START ===
    # Sortiere die Schlüssel, damit die Berechnungen der Reihe nach (Probe 1, Probe 2, Probe 3) ausgeführt werden
    for probe_name in sort(collect(keys(file_paths)))
        filepath = file_paths[probe_name]
        fit_and_plot_peaks(probe_name, filepath, json_peak_data)
    end
    # === GEÄNDERTE STELLE ENDE ===

    println("\nAlle Fitting- und Plotting-Prozesse abgeschlossen.")
    println("Prüfen Sie die Ordner 'Gefittete_Diffraktogramme_PDFs' und 'Fittergebnisse_CSV' für die Ergebnisse.")
end

# Führe den gesamten Fitting-Prozess aus
main_fitting_process()
