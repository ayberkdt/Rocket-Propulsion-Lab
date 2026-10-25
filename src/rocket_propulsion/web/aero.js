/* Rocket Propulsion Lab — gerçek akışkan ve aerodinamik modülleri.
 *
 * Bu dosya hesap yapmaz: her sayı /api/v1/real-fluid, /airfoil, /wing,
 * /supersonic-body, /flat-plate, /entry ve /launch-vehicle-aero uçlarından gelir
 * (src/rocket_propulsion/api/aero_routes.py). Burada yalnız birim gösterimi
 * (radyan → derece) yapılır ve grafikler çizilir. Genel işlevler app.js'ten
 * (calculate, renderResults, renderLineChart, formatValue) ve visuals.js'ten
 * (window.labCharts) gelir.                                                      */

(() => {
  const byId = (id) => document.getElementById(id);
  const degrees = (radians) => (radians * 180) / Math.PI;
  const percent = (value) => (value === null || value === undefined
    ? "çözülmedi" : `${(value * 100).toFixed(2).replace(".", ",")} %`);

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

  const renderers = {
    "real-fluid-form": ["/api/v1/real-fluid", (data) => {
      scoped("real-fluid", {
        ...data.state,
        cubic_difference: percent(data.cubic_relative_density_difference),
      });
      // The isobar jumps at the saturation temperature: draw each phase as its own
      // series so no line bridges the two-phase gap (no intermediate states exist).
      const liquid = data.isobar.filter((row) => row.phase === "liquid");
      const fluid = data.isobar.filter((row) => row.phase !== "liquid");
      line("real-fluid-density-chart", [
        [liquid, "temperature_k", "density_kg_m3", 1], [fluid, "temperature_k", "density_kg_m3", 2],
      ], { xLabel: "T (K)", yLabel: "ρ (kg/m³)" });
      line("real-fluid-cp-chart", [
        [liquid, "temperature_k", "cp_j_kg_k", 1], [fluid, "temperature_k", "cp_j_kg_k", 2],
      ], { xLabel: "T (K)", yLabel: "cₚ (J/(kg·K))" });
      line("real-fluid-saturation-chart", [[data.saturation, "temperature_k", "pressure_pa", 3]], { xLabel: "T (K)", yLabel: "p_sat (Pa)" });
      warn("real-fluid-warning", data.warnings || []);
    }],

    "airfoil-form": ["/api/v1/airfoil", (data) => {
      scoped("airfoil-solver", {
        panel_lift_coefficient: data.panel.lift_coefficient,
        thin_lift_coefficient: data.thin.lift_coefficient,
        zero_lift_angle_deg: degrees(data.thin.zero_lift_angle_rad),
        thin_cm_quarter: data.thin.moment_coefficient_quarter_chord,
        center_of_pressure: data.thin.center_of_pressure,
        panel_count: data.panel.panel_count,
      });
      // Aerodynamic convention: suction (negative Cp) plotted upward.
      const flipped = data.cp.map((row) => ({ ...row, minus_cp: -row.cp }));
      line("airfoil-cp-chart", [
        [flipped.filter((row) => row.upper), "x_over_c", "minus_cp", 1],
        [flipped.filter((row) => !row.upper), "x_over_c", "minus_cp", 2],
      ], { xLabel: "x/c", yLabel: "−Cp" });
      line("airfoil-lift-chart", [
        [data.sweep, "angle_of_attack_deg", "panel_lift_coefficient", 1, { points: true }],
        [data.sweep, "angle_of_attack_deg", "thin_lift_coefficient", "reference"],
      ], { xLabel: "α (°)", yLabel: "c_l" });
      warn("airfoil-warning", data.thin.warnings || []);
    }],

    "wing-form": ["/api/v1/wing", (data) => {
      scoped("wing-solver", data);
      line("wing-span-chart", [[data.spanwise, "y_over_semispan", "section_lift_coefficient", 1]], { xLabel: "2y/b", yLabel: "c_l" });
      line("wing-taper-chart", [[data.sweep, "taper_ratio", "induced_drag_factor", 2, { points: true }]], { xLabel: "λ", yLabel: "δ" });
    }],

    "supersonic-form": ["/api/v1/supersonic-body", (data) => {
      scoped("supersonic-solver", { ...data.cone, ...data.at_angle });
      line("supersonic-lift-chart", [
        [data.sweep, "angle_of_attack_deg", "shock_expansion_lift", 1],
        [data.sweep, "angle_of_attack_deg", "busemann_lift", 2],
        [data.sweep, "angle_of_attack_deg", "ackeret_lift", "reference"],
      ], { xLabel: "α (°)", yLabel: "c_l" });
      line("supersonic-drag-chart", [
        [data.sweep, "angle_of_attack_deg", "shock_expansion_drag", 1],
        [data.sweep, "angle_of_attack_deg", "busemann_drag", 2],
        [data.sweep, "angle_of_attack_deg", "ackeret_drag", "reference"],
      ], { xLabel: "α (°)", yLabel: "c_d" });
    }],

    "flat-plate-form": ["/api/v1/flat-plate", (data) => {
      scoped("flat-plate", {
        van_driest_cf: data.van_driest_ii.average_skin_friction,
        laminar_cf: data.laminar.average_skin_friction,
        turbulent_cf: data.turbulent.average_skin_friction,
        reynolds_length: data.reynolds_length,
        adiabatic_wall_temperature_k: data.turbulent.adiabatic_wall_temperature_k,
        edge_temperature_k: data.atmosphere.temperature_k,
      });
      line("flat-plate-cf-chart", [
        [data.sweep, "x_m", "turbulent_cf", 1],
        [data.sweep, "x_m", "laminar_cf", 2],
      ], { xLabel: "x (m)", yLabel: "c_f" });
      warn("flat-plate-warning", [...(data.laminar.warnings || []), ...(data.turbulent.warnings || [])]);
    }],

    "entry-form": ["/api/v1/entry", (data) => {
      scoped("entry-solver", {
        peak_deceleration_g: data.peak_deceleration.deceleration_g,
        peak_deceleration_altitude_m: data.peak_deceleration.altitude_m,
        peak_heating_w_m2: data.peak_heating ? data.peak_heating.stagnation_heating_w_m2 : null,
        integrated_heat_load_j_m2: data.integrated_heat_load_j_m2,
        termination: data.termination === "velocity" ? "hız eşiği (100 m/s)" : data.termination === "altitude" ? "yer" : "süre sınırı",
      });
      line("entry-altitude-chart", [[data.trajectory, "velocity_m_s", "altitude_m", 1]], { xLabel: "V (m/s)", yLabel: "h (m)" });
      line("entry-decel-chart", [[data.trajectory, "time_s", "deceleration_g", 4]], { xLabel: "t (s)", yLabel: "a (g)" });
      line("entry-heating-chart", [[data.trajectory, "time_s", "stagnation_heating_w_m2", 2]], { xLabel: "t (s)", yLabel: "q̇ (W/m²)" });
      warn("entry-warning", data.warnings || []);
    }],

    "launch-aero-form": ["/api/v1/launch-vehicle-aero", (data) => {
      scoped("launch-aero-solver", {
        ...data.normal_force,
        normal_force_coefficient: data.crossflow.normal_force_coefficient,
        dynamic_pressure_pa: data.dynamic_pressure_pa,
        normal_force_n: data.normal_force_n,
        friction_drag_coefficient: data.drag.friction_drag_coefficient,
        base_drag_coefficient: data.drag.base_drag_coefficient,
      });
      line("launch-cn-chart", [
        [data.sweep, "angle_of_attack_deg", "normal_force_coefficient", 1],
        [data.sweep, "angle_of_attack_deg", "potential_normal_force_coefficient", "reference"],
      ], { xLabel: "α (°)", yLabel: "C_N" });
      const drag = [
        { label: "Sürtünme", value: data.drag.friction_drag_coefficient, slot: 1 },
        { label: "Taban", value: data.drag.base_drag_coefficient, slot: 2 },
      ];
      if (data.drag.forebody_pressure_drag_coefficient !== null) {
        drag.push({ label: "Ön gövde basınç", value: data.drag.forebody_pressure_drag_coefficient, slot: 3 });
      }
      window.labCharts.renderBars(byId("launch-drag-chart"), drag.map((row) => ({
        label: row.label, segments: [{ value: row.value, slot: row.slot, name: "C_D" }],
      })), { unit: "", xLabel: "Sürükleme katsayısı C_D" });
      warn("launch-aero-warning", [...(data.crossflow.warnings || []), ...(data.drag.base.warnings || [])]);
    }],
  };

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
