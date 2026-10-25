/* Rocket Propulsion Lab — ısı transferi, rejeneratif soğutma ve yanma ürünleri.
 *
 * Bu dosya hesap yapmaz: her sayı /api/v1/heat-*, /regen-chamber, /transpiration-cooling,
 * /combustion-stability, /propellant-vaporization,
 * /combustion-products ve /combustion-losses uçlarından gelir
 * (api/thermal_routes.py, api/combustion_dynamics_routes.py). Burada yalnız birim gösterimi
 * (m → mm) ve çizim yapılır. Genel işlevler app.js'ten (calculate,
 * renderResults, renderLineChart) ve visuals.js'ten (window.labCharts) gelir.   */

(() => {
  const byId = (id) => document.getElementById(id);

  function scoped(moduleName, values) {
    renderResults(document.querySelector(`[data-module="${moduleName}"]`), values);
  }

  function line(id, series, labels) {
    renderLineChart(byId(id), series.map(([data, x, y, slot, extra = {}]) => ({
      data, x, y, className: `series-${slot}`, ...extra,
    })), labels);
  }

  function warn(id, messages) {
    const node = byId(id);
    if (node) node.textContent = messages.filter(Boolean).map((item) => item.message || item).join(" ");
  }

  function bars(id, rows, unit, xLabel) {
    window.labCharts.renderBars(byId(id), rows.map((row) => ({
      label: row.label, segments: [{ value: row.value, slot: row.slot || 1, name: row.name || "" }],
    })), { unit, xLabel });
  }

  const renderers = {
    "heat-conduction-form": ["/api/v1/heat-conduction", (data) => {
      scoped("heat-conduction", data);
      line("conduction-history-chart", [
        [data.history, "time_s", "hot_surface_k", 1], [data.history, "time_s", "cold_surface_k", 2],
      ], { xLabel: "t (s)", yLabel: "T (K)" });
      const picks = [data.profiles[0], data.profiles[Math.floor(data.profiles.length / 2)], data.profiles[data.profiles.length - 1]];
      line("conduction-profile-chart", picks.map((profile, index) => [profile.points, "x_mm", "temperature_k", [1, 3, 2][index]]),
        { xLabel: "x (mm, soğutucu yüzünden)", yLabel: "T (K)" });
      warn("conduction-warning", data.warnings || []);
    }],

    "heat-convection-form": ["/api/v1/heat-convection", (data) => {
      scoped("heat-convection", data.transport);
      line("convection-h-chart", [
        [data.sweep, "mass_flux_kg_m2_s", "gnielinski", 1],
        [data.sweep, "mass_flux_kg_m2_s", "dittus_boelter", 2],
        [data.sweep, "mass_flux_kg_m2_s", "colburn", "reference"],
      ], { xLabel: "G (kg/(m²·s))", yLabel: "h (W/(m²·K))" });
      warn("convection-warning", data.transport.warnings || []);
    }],

    "heat-radiation-form": ["/api/v1/heat-radiation", (data) => {
      scoped("heat-radiation", data);
      line("radiation-emissivity-chart", [[data.sweep, "emissivity", "wall_temperature_k", 1, { points: true }]],
        { xLabel: "ε", yLabel: "T_w (K)" });
    }],

    "regen-form": ["/api/v1/regen-chamber", (data) => {
      scoped("regen-solver", data);
      const stations = data.stations.map((row) => ({ ...row, x_mm: row.x_m * 1000 }));
      line("regen-temperature-chart", [
        [stations, "x_mm", "hot_wall_temperature_k", 1],
        [stations, "x_mm", "cold_wall_temperature_k", 2],
        [stations, "x_mm", "coolant_temperature_k", 3],
      ], { xLabel: "x (mm, boğaz = 0)", yLabel: "T (K)" });
      line("regen-flux-chart", [[stations, "x_mm", "heat_flux_w_m2", 4]], { xLabel: "x (mm)", yLabel: "q (W/m²)" });
      line("regen-pressure-chart", [[stations, "x_mm", "coolant_pressure_pa", 5]], { xLabel: "x (mm)", yLabel: "p (Pa)" });
      warn("regen-warning", data.warnings || []);
    }],

    "combustion-products-form": ["/api/v1/combustion-products", (data) => {
      scoped("combustion-products", {
        ...data.dissociation,
        vacuum_specific_impulse_s: data.vacuum_specific_impulse_s,
        stoichiometric_oxidizer_fuel_ratio: data.stoichiometry.stoichiometric_oxidizer_fuel_ratio,
      });
      const chamber = data.stations[0];
      bars("products-chamber-chart", chamber.leading_species.slice(0, 8).map((row) => ({
        label: row.species.replace(/\d/g, (digit) => "₀₁₂₃₄₅₆₇₈₉"[digit]), value: row.mole_fraction, slot: 1, name: "mol kesri",
      })), "", "Mol kesri (oda)");
      const names = { chamber: "Oda", throat: "Boğaz", exit: "Çıkış" };
      bars("products-radical-chart", data.stations.map((station, index) => ({
        label: names[station.name] || station.name, value: station.radical_mole_fraction, slot: [4, 2, 3][index] || 1, name: "radikal",
      })), "", "Toplam radikal mol kesri");
      warn("products-warning", data.warnings || []);
    }],

    "combustion-losses-form": ["/api/v1/combustion-losses", (data) => {
      scoped("combustion-losses", { ...data.chain, ...data.finite_area_combustor });
      const chain = data.chain;
      bars("losses-chain-chart", [
        { label: "Kinetik (ODK)", value: chain.kinetic_loss_s, slot: 1 },
        { label: "Iraksama (2B)", value: chain.divergence_loss_s, slot: 2 },
        { label: "Sınır tabaka", value: chain.boundary_layer_loss_s, slot: 3 },
        { label: "Enerji salımı (c*)", value: chain.energy_release_loss_s, slot: 4 },
      ].map((row) => ({ ...row, name: "kayıp" })), "s", "Isp kaybı (s)");
      const kinetics = data.kinetics ? data.kinetics.stations : [];
      line("losses-radicals-chart", kinetics.length ? [
        [kinetics, "area_ratio", "x_H", 1], [kinetics, "area_ratio", "x_OH", 2],
      ] : [], { xLabel: "A/A*", yLabel: "mol kesri" });
      warn("losses-warning", data.kinetics ? data.kinetics.warnings : ["Sonlu hızlı kimya yalnız LOX/LH₂ için; kinetik verim 1 alındı."]);
    }],
  };

  renderers["transpiration-form"] = ["/api/v1/transpiration-cooling", (data) => {
    scoped("transpiration", data.required);
    line("transpiration-wall-chart", [[data.sweep, "coolant_mass_flux_kg_m2_s", "wall_temperature_k", 1]],
      { xLabel: "ṁ″ (kg/(m²·s))", yLabel: "duvar T (K)" });
    line("transpiration-blowing-chart", [[data.sweep, "coolant_mass_flux_kg_m2_s", "blowing_correction", 2]],
      { xLabel: "ṁ″ (kg/(m²·s))", yLabel: "h/h*" });
    warn("transpiration-warning", data.warnings || []);
  }];

  renderers["combustion-stability-form"] = ["/api/v1/combustion-stability", (data) => {
    const chug = data.chugging;
    const oneT = data.modes.find((mode) => mode.name === "1T");
    scoped("combustion-stability", {
      verdict: chug.unconditionally_stable ? "koşulsuz kararlı (K ≤ 1)" : (chug.stable ? "kararlı" : "kararsız"),
      gain: chug.gain,
      critical_time_lag_s: chug.critical_time_lag_s,
      neutral_frequency_hz: chug.neutral_frequency_hz,
      residence_time_s: data.chamber.residence_time_s,
      first_tangential_hz: oneT ? oneT.frequency_hz : null,
    });
    bars("stability-modes-chart", data.modes.slice(0, 8).map((mode, index) => ({
      label: mode.name, value: mode.frequency_hz, slot: mode.name.includes("T") || mode.name.includes("R") ? 2 : 1,
      name: mode.name.includes("T") || mode.name.includes("R") ? "enine/karma" : "boyuna",
    })), "Hz", "Frekans (Hz)");
    line("stability-nyquist-chart", [[data.nyquist.curve, "real", "imaginary", 1]],
      { xLabel: "Re f(iω)", yLabel: "Im f(iω)" });
    const zeros = data.nyquist.right_half_plane_zeros;
    warn("stability-warning", [`Nyquist: sağ yarı düzlemde ${zeros} kök${zeros ? " (kararsız)" : ""}; −1 noktası çevrelenmemeli. Rijit besleme; besleme hattı rezonansları için Python API'deki feed_coupled_stability kullanılır.`]);
  }];

  renderers["propellant-vaporization-form"] = ["/api/v1/propellant-vaporization", (data) => {
    scoped("propellant-vaporization", { ...data.result, oxidizer_generalized_length_m: data.oxidizer_generalized_length_m,
      fuel_generalized_length_m: data.fuel_generalized_length_m });
    line("vaporization-length-chart", [
      [data.sweep, "chamber_length_m", "oxidizer_percent_vaporized", 1],
      [data.sweep, "chamber_length_m", "fuel_percent_vaporized", 2],
    ], { xLabel: "L (m)", yLabel: "buharlaşan kütle (%)" });
    warn("vaporization-warning", data.properties.fuel ? [] : ["Hidrojen gaz olarak enjekte edilir; tamamen buharlaşmış sayılır."]);
  }];

  for (const [formId, [endpoint, render]] of Object.entries(renderers)) {
    const form = byId(formId);
    if (!form) continue;
    form.addEventListener("submit", (event) => {
      event.preventDefault();
      calculate(form, endpoint, render);
    });
    form.requestSubmit();
  }
})();
