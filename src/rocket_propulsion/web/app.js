const numberFormatter = new Intl.NumberFormat("tr-TR", {
  maximumSignificantDigits: 9,
  useGrouping: false,
});

function formatValue(value, key) {
  if (value === null || value === undefined) return "tanımsız";
  if (typeof value === "string") {
    const labels = {
      adapted: "uyumlu",
      underexpanded: "az genişlemiş",
      overexpanded: "aşırı genişlemiş",
      compression: "sıkıştırma",
      expansion: "genleşme",
      "no-flow": "akış yok",
      unchoked: "boğulmamış",
      "just-choked": "yeni boğulmuş",
      "internal-normal-shock": "nozul içinde normal şok",
      "shock-at-exit": "çıkışta normal şok",
      design: "tasarım noktası",
      target_reached: "hedefe ulaşıldı",
      tank_reserve: "tank rezervi",
      dry_mass_floor: "kuru kütle sınırı",
      maximum_duration: "maksimum süre",
    };
    return labels[value] || value;
  }
  if (key.endsWith("_deg")) return `${numberFormatter.format(value)}°`;
  const magnitude = Math.abs(value);
  if ((magnitude > 0 && magnitude < 1e-5) || magnitude >= 1e7) {
    return Number(value).toExponential(6);
  }
  return numberFormatter.format(value);
}

function formPayload(form) {
  return Object.fromEntries(new FormData(form).entries());
}

function apiErrorMessage(payload, fallback) {
  if (typeof payload?.error === "string") return payload.error;
  if (payload?.error?.message) return payload.error.message;
  return fallback;
}

function renderResults(module, data) {
  module.querySelectorAll("[data-key]").forEach((element) => {
    const key = element.dataset.key;
    element.textContent = formatValue(data[key], key);
  });
}

const SVG_NAMESPACE = "http://www.w3.org/2000/svg";

function svgElement(tag, attributes = {}, text = "") {
  const element = document.createElementNS(SVG_NAMESPACE, tag);
  for (const [name, value] of Object.entries(attributes)) element.setAttribute(name, value);
  if (text) element.textContent = text;
  return element;
}

function chartTick(value) {
  const magnitude = Math.abs(value);
  if (magnitude >= 10000 || (magnitude > 0 && magnitude < 0.001)) {
    return value.toExponential(1);
  }
  return Number(value.toPrecision(4)).toString();
}

function renderLineChart(svg, series, { xLabel, yLabel }) {
  const width = 720;
  const height = 280;
  const margin = { top: 16, right: 20, bottom: 48, left: 72 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const allPoints = series.flatMap((item) => item.data.map((point) => ({
    x: Number(point[item.x]),
    y: Number(point[item.y]),
  }))).filter((point) => Number.isFinite(point.x) && Number.isFinite(point.y));
  svg.replaceChildren();
  if (!allPoints.length) return;

  let xMin = Math.min(...allPoints.map((point) => point.x));
  let xMax = Math.max(...allPoints.map((point) => point.x));
  let yMin = Math.min(...allPoints.map((point) => point.y));
  let yMax = Math.max(...allPoints.map((point) => point.y));
  const xPadding = Math.max((xMax - xMin) * 0.03, 1e-12);
  const yPadding = Math.max((yMax - yMin) * 0.08, 1e-12);
  xMin -= xPadding;
  xMax += xPadding;
  yMin -= yPadding;
  yMax += yPadding;

  const xScale = (value) => margin.left + ((value - xMin) / (xMax - xMin)) * plotWidth;
  const yScale = (value) => margin.top + plotHeight - ((value - yMin) / (yMax - yMin)) * plotHeight;
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);

  for (let index = 0; index <= 5; index += 1) {
    const fraction = index / 5;
    const x = margin.left + fraction * plotWidth;
    const y = margin.top + fraction * plotHeight;
    svg.append(
      svgElement("line", { x1: x, y1: margin.top, x2: x, y2: margin.top + plotHeight, class: "chart-grid-line" }),
      svgElement("line", { x1: margin.left, y1: y, x2: margin.left + plotWidth, y2: y, class: "chart-grid-line" }),
      svgElement("text", { x, y: margin.top + plotHeight + 20, "text-anchor": "middle", class: "chart-tick" }, chartTick(xMin + fraction * (xMax - xMin))),
      svgElement("text", { x: margin.left - 10, y: y + 4, "text-anchor": "end", class: "chart-tick" }, chartTick(yMax - fraction * (yMax - yMin))),
    );
  }
  svg.append(
    svgElement("line", { x1: margin.left, y1: margin.top + plotHeight, x2: margin.left + plotWidth, y2: margin.top + plotHeight, class: "chart-axis" }),
    svgElement("line", { x1: margin.left, y1: margin.top, x2: margin.left, y2: margin.top + plotHeight, class: "chart-axis" }),
    svgElement("text", { x: margin.left + plotWidth / 2, y: height - 8, "text-anchor": "middle", class: "chart-axis-title" }, xLabel),
    svgElement("text", { x: 16, y: margin.top + plotHeight / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${margin.top + plotHeight / 2})`, class: "chart-axis-title" }, yLabel),
  );

  series.forEach((item, index) => {
    const path = item.data
      .map((point, pointIndex) => `${pointIndex ? "L" : "M"}${xScale(Number(point[item.x])).toFixed(2)},${yScale(Number(point[item.y])).toFixed(2)}`)
      .join(" ");
    svg.append(svgElement("path", { d: path, class: `chart-line ${item.className || `series-${index + 1}`}` }));
  });
}

function renderPolytropicCharts(data) {
  renderLineChart(
    document.getElementById("pv-chart"),
    [{ data: data.path, x: "specific_volume_m3_kg", y: "pressure_pa", className: "series-1" }],
    { xLabel: "v (m³/kg)", yLabel: "p (Pa)" },
  );
  renderLineChart(
    document.getElementById("ts-chart"),
    [{ data: data.path, x: "entropy_change_j_kg_k", y: "temperature_k", className: "series-2" }],
    { xLabel: "Δs (J/kg·K)", yLabel: "T (K)" },
  );
}

let latestGeometryExports = null;

function renderGeometryComparison(records = []) {
  const body = document.getElementById("geometry-comparison-body");
  body.replaceChildren();
  records.forEach((record) => {
    const row = document.createElement("tr");
    [
      record.contour,
      formatValue(record.diverging_length_m, "diverging_length_m"),
      formatValue(record.surface_area_m2, "surface_area_m2"),
      formatValue(record.estimated_divergence_efficiency, "estimated_divergence_efficiency"),
      formatValue(record.exit_flow_angle_deg, "exit_flow_angle_deg"),
    ].forEach((value) => row.append(textElement("td", value)));
    body.append(row);
  });
}

function renderNozzleGeometry(data) {
  const mirrored = data.stations.map((point) => ({ ...point, mirrored_radius_m: -point.radius_m }));
  renderLineChart(
    document.getElementById("nozzle-characteristic-chart"),
    (data.characteristic_lines || []).map((line) => ({
      data: line,
      x: "x_m",
      y: "radius_m",
      className: "series-1",
    })),
    { xLabel: "x (m)", yLabel: "r (m)" },
  );
  renderLineChart(
    document.getElementById("nozzle-contour-chart"),
    [
      { data: data.stations, x: "x_m", y: "radius_m", className: "series-neutral" },
      { data: mirrored, x: "x_m", y: "mirrored_radius_m", className: "series-neutral" },
    ],
    { xLabel: "x (m)", yLabel: "r (m)" },
  );
  renderLineChart(
    document.getElementById("nozzle-mach-chart"),
    [{ data: data.stations, x: "x_m", y: "mach", className: "series-1" }],
    { xLabel: "x (m)", yLabel: "Mach" },
  );
  renderLineChart(
    document.getElementById("nozzle-state-chart"),
    [
      { data: data.stations, x: "x_m", y: "pressure_total_ratio", className: "series-1" },
      { data: data.stations, x: "x_m", y: "temperature_total_ratio", className: "series-2" },
    ],
    { xLabel: "x (m)", yLabel: "Statik / toplam" },
  );
  latestGeometryExports = data.exports || null;
  renderGeometryComparison(data.comparison);
}

async function downloadGeometryExport(format) {
  if (!latestGeometryExports?.[format]) {
    const form = document.getElementById("geometry-form");
    const payload = formPayload(form);
    payload.include_exports = "true";
    payload.include_comparison = "false";
    const response = await fetch("/api/v1/nozzle-contour", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const result = await response.json();
    if (!response.ok) throw new Error(apiErrorMessage(result, "Dışa aktarma üretilemedi."));
    latestGeometryExports = result.data.exports;
  }
  const mediaTypes = {
    csv: "text/csv;charset=utf-8",
    svg: "image/svg+xml;charset=utf-8",
    dxf: "application/dxf;charset=utf-8",
    stl: "model/stl;charset=utf-8",
  };
  const blob = new Blob([latestGeometryExports[format]], { type: mediaTypes[format] });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = `rocket-nozzle.${format}`;
  link.click();
  URL.revokeObjectURL(link.href);
}

function renderLossBudget(data) {
  document.getElementById("loss-band-low").textContent = formatValue(data.delivered_thrust_band_n[0], "thrust_n");
  document.getElementById("loss-band-high").textContent = formatValue(data.delivered_thrust_band_n[1], "thrust_n");
  const terms = data.terms.map((term, index) => ({ ...term, index: index + 1 }));
  renderLineChart(
    document.getElementById("loss-budget-chart"),
    [{ data: terms, x: "index", y: "thrust_loss_n", className: "series-2" }],
    { xLabel: "Terim sırası", yLabel: "İtki kaybı (N)" },
  );
}

function renderTwoPhaseBand(data) {
  document.getElementById("two-phase-thrust-low").textContent = formatValue(data.thrust_range_n[0], "thrust_n");
  document.getElementById("two-phase-thrust-high").textContent = formatValue(data.thrust_range_n[1], "thrust_n");
  document.getElementById("two-phase-isp-low").textContent = formatValue(data.specific_impulse_range_s[0], "specific_impulse_s");
  document.getElementById("two-phase-isp-high").textContent = formatValue(data.specific_impulse_range_s[1], "specific_impulse_s");
}

function renderNozzleThermal(data) {
  const cooling = data.cooling || {};
  document.getElementById("coolant-outlet-temperature").textContent = formatValue(cooling.outlet_temperature_k, "temperature_k");
  document.getElementById("coolant-pressure-drop").textContent = formatValue(cooling.pressure_drop_pa, "pressure_pa");
  document.getElementById("coolant-energy-residual").textContent = formatValue(cooling.energy_balance_residual_w, "residual_w");
  renderLineChart(
    document.getElementById("heat-flux-chart"),
    [{ data: data.stations, x: "x_m", y: "heat_flux_w_m2", className: "series-2" }],
    { xLabel: "x (m)", yLabel: "q'' (W/m²)" },
  );
  renderLineChart(
    document.getElementById("adiabatic-wall-chart"),
    [{ data: data.stations, x: "x_m", y: "adiabatic_wall_temperature_k", className: "series-1" }],
    { xLabel: "x (m)", yLabel: "Taw (K)" },
  );
  const warnings = [...(data.warnings || []), ...(cooling.warnings || [])];
  document.getElementById("thermal-warning").textContent = warnings.length
    ? warnings.map((warning) => warning.message).join(" ")
    : "Bartz ve kanal enerji dengesi tanımlı girdi aralıklarında kapandı.";
}

let latestStagingResult = null;
let latestCycleResult = null;
let latestUncertaintyResult = null;
let latestBurnInput = null;
let latestBurnResult = null;
let loadedWorkspace = null;

function burnRequest(form) {
  const values = formPayload(form);
  const request = {
    provider: {
      family: "bipropellant",
      name: "UI constant bipropellant point",
      delivered_thrust_n: values.delivered_thrust_n,
      ideal_thrust_n: values.ideal_thrust_n,
      system_specific_impulse_s: values.system_specific_impulse_s,
      chamber_specific_impulse_s: values.chamber_specific_impulse_s,
      oxidizer_fuel_ratio: values.oxidizer_fuel_ratio,
      dump_mass_flow_kg_s: values.dump_mass_flow_kg_s,
      dump_role: values.dump_role,
      chamber_pressure_pa: values.chamber_pressure_pa,
      ambient_pressure_pa: 0,
      oxidizer_tank_id: "oxidizer",
      fuel_tank_id: "fuel",
      oxidizer_key: "oxidizer",
      fuel_key: "fuel",
    },
    definition: {
      name: values.profile_mode === "piecewise" ? "UI profiled burn" : "UI constant burn",
      initial_mass_kg: values.initial_mass_kg,
      protected_dry_mass_kg: values.protected_dry_mass_kg,
      maximum_duration_s: values.maximum_duration_s,
      cant_efficiency: values.cant_efficiency,
      target: { kind: values.target_kind, value: values.target_value },
      tanks: [
        {
          tank_id: "oxidizer",
          propellant_key: "oxidizer",
          role: "oxidizer",
          loaded_mass_kg: values.oxidizer_loaded_mass_kg,
          reserve_mass_kg: values.oxidizer_reserve_mass_kg,
        },
        {
          tank_id: "fuel",
          propellant_key: "fuel",
          role: "fuel",
          loaded_mass_kg: values.fuel_loaded_mass_kg,
          reserve_mass_kg: values.fuel_reserve_mass_kg,
        },
      ],
    },
  };
  if (values.profile_mode === "piecewise") {
    const segments = [];
    if (values.ramp_up_s > 0) {
      segments.push({ duration_s: values.ramp_up_s, start_throttle: 0, end_throttle: 1, phase: "ramp-up" });
    }
    if (values.steady_duration_s > 0) {
      segments.push({ duration_s: values.steady_duration_s, start_throttle: 1, end_throttle: 1, phase: "steady" });
    }
    if (values.ramp_down_s > 0) {
      segments.push({ duration_s: values.ramp_down_s, start_throttle: 1, end_throttle: 0, phase: "ramp-down" });
    }
    request.schedule = { name: "UI ramp-steady-ramp schedule", segments };
  }
  return request;
}

function renderBurn(data, request = null) {
  latestBurnResult = data;
  if (request) latestBurnInput = request;
  const module = document.querySelector('[data-module="burn"]');
  renderResults(module, {
    ...data.summary,
    fidelity_level: data.schema === "rocket_propulsion_burn_v2" ? "L1 · profilli" : "L0 · sabit",
  });
  renderLineChart(
    document.getElementById("burn-mass-chart"),
    [{ data: data.profile, x: "time_s", y: "vehicle_mass_kg", className: "series-1" }],
    { xLabel: "t (s)", yLabel: "Araç kütlesi (kg)" },
  );
  renderLineChart(
    document.getElementById("burn-thrust-chart"),
    [{ data: data.profile, x: "time_s", y: "delivered_thrust_n", className: "series-2" }],
    { xLabel: "t (s)", yLabel: "İtki (N)" },
  );
  const finalTanks = Object.fromEntries(data.final_state.tank_masses_kg);
  const table = document.getElementById("burn-tank-body");
  table.replaceChildren();
  data.summary.consumed_by_tank_kg.forEach(([tankId, consumed]) => {
    const row = document.createElement("tr");
    [tankId, formatValue(consumed, "mass_kg"), formatValue(finalTanks[tankId], "mass_kg")]
      .forEach((value) => row.append(textElement("td", value)));
    table.append(row);
  });
  ["burn-json", "burn-report", "burn-sidera"].forEach((id) => {
    document.getElementById(id).disabled = false;
  });
  const workspaceMessage = document.getElementById("workspace-message");
  if (workspaceMessage) workspaceMessage.textContent = "Burn vakası çalışma dosyasına eklenmeye hazır.";
}

async function runBurn(form) {
  const message = form.querySelector(".form-message");
  const button = form.querySelector("button[type='submit']");
  message.textContent = "Motor çalışma noktası çözülüyor…";
  message.className = "form-message";
  button.disabled = true;
  try {
    const request = burnRequest(form);
    const pointResponse = await fetch("/api/v1/engine-operating-point", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(request.provider),
    });
    const pointPayload = await pointResponse.json();
    if (!pointResponse.ok) throw new Error(apiErrorMessage(pointPayload, "Çalışma noktası çözülemedi."));
    message.textContent = "Tank envanteri ve burn hedefi çözülüyor…";
    const simulationRequest = {
      definition: request.definition,
      operating_point: pointPayload.data,
    };
    if (request.schedule) simulationRequest.schedule = request.schedule;
    const response = await fetch("/api/v1/burn/simulate", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(simulationRequest),
    });
    const payload = await response.json();
    if (!response.ok) throw new Error(apiErrorMessage(payload, "Burn çözümü tamamlanamadı."));
    renderBurn(payload.data, simulationRequest);
    message.textContent = payload.data.summary.target_achieved
      ? "Burn hedefi analitik olarak kapandı."
      : `Hedefe ulaşılamadı: ${formatValue(payload.data.summary.cutoff_reason, "cutoff_reason")}.`;
    message.className = payload.data.summary.target_achieved
      ? "form-message success"
      : "form-message error";
  } catch (error) {
    message.textContent = error.message;
    message.className = "form-message error";
  } finally {
    button.disabled = false;
  }
}

function renderStaging(data) {
  latestStagingResult = data;
  document.getElementById("workspace-message").textContent = "Sistem analizleri çalışma dosyasına eklenmeye hazır.";
  const stages = data.stages.map((stage, index) => ({ ...stage, index: index + 1 }));
  renderLineChart(
    document.getElementById("staging-chart"),
    [{ data: stages, x: "index", y: "delta_v_m_s", className: "series-3" }],
    { xLabel: "Kademe sırası", yLabel: "Δv (m/s)" },
  );
}

function renderEngineCycle(data) {
  latestCycleResult = data;
  document.getElementById("workspace-message").textContent = "Sistem analizleri çalışma dosyasına eklenmeye hazır.";
}

function renderUncertainty(data) {
  latestUncertaintyResult = data;
  document.getElementById("workspace-message").textContent = "Sistem analizleri çalışma dosyasına eklenmeye hazır.";
  document.getElementById("uncertainty-p05").textContent = formatValue(data.output_percentiles.p05, "delta_v_m_s");
  document.getElementById("uncertainty-p95").textContent = formatValue(data.output_percentiles.p95, "delta_v_m_s");
  const samples = data.samples.map((sample, index) => ({ index: index + 1, output: sample.output }));
  renderLineChart(
    document.getElementById("uncertainty-chart"),
    [{ data: samples, x: "index", y: "output", className: "series-1" }],
    { xLabel: "Örnek", yLabel: "Δv (m/s)" },
  );
}

function analysisWorkspace() {
  const cases = [];
  if (latestStagingResult) cases.push({
    case_id: "staging",
    name: "Kademe optimizasyonu",
    model: "staging-optimizer",
    inputs: formPayload(document.getElementById("staging-form")),
    results: {
      initial_mass_kg: latestStagingResult.initial_mass_kg,
      payload_fraction: latestStagingResult.payload_fraction,
      achieved_delta_v_m_s: latestStagingResult.achieved_delta_v_m_s,
    },
    model_versions: { "staging-optimizer": "1.0" },
    warnings: latestStagingResult.active_constraints || [],
  });
  if (latestCycleResult) cases.push({
    case_id: "engine-cycle",
    name: "Motor çevrimi",
    model: latestCycleResult.cycle,
    inputs: formPayload(document.getElementById("engine-cycle-form")),
    results: {
      pump_power_w: latestCycleResult.pump_power_w,
      turbine_power_w: latestCycleResult.turbine_power_w,
      delivered_mass_fraction: latestCycleResult.delivered_mass_fraction,
    },
    model_versions: { "engine-cycle": "1.0" },
    warnings: (latestCycleResult.warnings || []).map((warning) => warning.message),
  });
  if (latestUncertaintyResult) cases.push({
    case_id: "uncertainty",
    name: "Δv belirsizliği",
    model: latestUncertaintyResult.method,
    inputs: formPayload(document.getElementById("uncertainty-form")),
    results: {
      output_mean: latestUncertaintyResult.output_mean,
      output_standard_deviation: latestUncertaintyResult.output_standard_deviation,
      p05: latestUncertaintyResult.output_percentiles.p05,
      p95: latestUncertaintyResult.output_percentiles.p95,
    },
    model_versions: { "uncertainty": "1.0" },
    warnings: [],
  });
  if (latestBurnResult && latestBurnInput) cases.push({
    case_id: "standalone-burn",
    name: latestBurnResult.definition.name,
    model: "standalone-constant-burn",
    inputs: latestBurnInput,
    results: {
      artifact: latestBurnResult,
      duration_s: latestBurnResult.summary.duration_s,
      consumed_propellant_kg: latestBurnResult.summary.consumed_propellant_kg,
      delivered_total_impulse_n_s: latestBurnResult.summary.delivered_total_impulse_n_s,
      ideal_delta_v_m_s: latestBurnResult.summary.ideal_delta_v_m_s,
      result_hash: latestBurnResult.result_hash,
    },
    model_versions: { "standalone-constant-burn": "analytic_constant_v1" },
    warnings: latestBurnResult.warnings || [],
  });
  return {
    schema_version: 2,
    name: "Rocket Propulsion Lab çalışma alanı",
    cases,
    dataset_versions: { thermochemistry: "NASA-Glenn-2002" },
    notes: "Tarayıcıdan oluşturuldu.",
  };
}

function downloadText(filename, content, mediaType) {
  const link = document.createElement("a");
  link.href = URL.createObjectURL(new Blob([content], { type: mediaType }));
  link.download = filename;
  link.click();
  setTimeout(() => URL.revokeObjectURL(link.href), 0);
}

async function normalizedWorkspace(workspace, includeReport = false) {
  const response = await fetch("/api/v1/workspace", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ workspace_json: JSON.stringify(workspace), include_report: includeReport }),
  });
  const payload = await response.json();
  if (!response.ok) throw new Error(apiErrorMessage(payload, "Çalışma dosyası doğrulanamadı."));
  return payload.data;
}

function renderThermalProperties(data) {
  const module = document.querySelector('[data-module="thermal-properties"]');
  renderResults(module, data.properties);
  renderLineChart(
    document.getElementById("thermal-cp-chart"),
    [
      { data: data.curve, x: "temperature_k", y: "cp_j_kg_k", className: "series-1" },
      { data: data.curve, x: "temperature_k", y: "cv_j_kg_k", className: "series-2" },
      { data: data.calorically_perfect_reference.curve, x: "temperature_k", y: "cp_j_kg_k", className: "series-reference" },
    ],
    { xLabel: "T (K)", yLabel: "J/(kg·K)" },
  );
  renderLineChart(
    document.getElementById("thermal-gamma-chart"),
    [
      { data: data.curve, x: "temperature_k", y: "gamma", className: "series-1" },
      { data: data.calorically_perfect_reference.curve, x: "temperature_k", y: "gamma", className: "series-reference" },
    ],
    { xLabel: "T (K)", yLabel: "γ" },
  );
  renderLineChart(
    document.getElementById("thermal-enthalpy-chart"),
    [
      { data: data.curve, x: "temperature_k", y: "enthalpy_j_kg", className: "series-3" },
      { data: data.calorically_perfect_reference.curve, x: "temperature_k", y: "enthalpy_j_kg", className: "series-reference" },
    ],
    { xLabel: "T (K)", yLabel: "h (J/kg)" },
  );
}

function renderEquilibriumSweep(data) {
  const module = document.querySelector('[data-module="equilibrium-sweep"]');
  const points = data.points.map((point) => ({
    ...point,
    molecular_mass_g_mol: point.molecular_mass_kg_mol * 1000,
  }));
  const peak = points.reduce((best, point) => (
    point.equilibrium_vacuum_specific_impulse_s > best.equilibrium_vacuum_specific_impulse_s
      ? point
      : best
  ));
  renderResults(module, {
    peak_ratio: peak.oxidizer_fuel_ratio,
    peak_temperature_k: peak.adiabatic_flame_temperature_k,
    peak_cstar_m_s: peak.characteristic_velocity_m_s,
    peak_frozen_isp_s: peak.frozen_vacuum_specific_impulse_s,
    peak_equilibrium_isp_s: peak.equilibrium_vacuum_specific_impulse_s,
    peak_element_residual: peak.maximum_element_relative_residual,
    provider_label: `${data.provider} · ${data.provider_version}`,
  });
  renderLineChart(
    document.getElementById("of-temperature-chart"),
    [{ data: points, x: "oxidizer_fuel_ratio", y: "adiabatic_flame_temperature_k", className: "series-1" }],
    { xLabel: "O/F", yLabel: "T₍c₎ (K)" },
  );
  renderLineChart(
    document.getElementById("of-cstar-chart"),
    [{ data: points, x: "oxidizer_fuel_ratio", y: "characteristic_velocity_m_s", className: "series-2" }],
    { xLabel: "O/F", yLabel: "c* (m/s)" },
  );
  renderLineChart(
    document.getElementById("of-molecular-mass-chart"),
    [{ data: points, x: "oxidizer_fuel_ratio", y: "molecular_mass_g_mol", className: "series-1" }],
    { xLabel: "O/F", yLabel: "M (g/mol)" },
  );
  renderLineChart(
    document.getElementById("of-gamma-chart"),
    [{ data: points, x: "oxidizer_fuel_ratio", y: "gamma", className: "series-3" }],
    { xLabel: "O/F", yLabel: "γ" },
  );
  renderLineChart(
    document.getElementById("of-isp-chart"),
    [
      { data: points, x: "oxidizer_fuel_ratio", y: "equilibrium_vacuum_specific_impulse_s", className: "series-1" },
      { data: points, x: "oxidizer_fuel_ratio", y: "frozen_vacuum_specific_impulse_s", className: "series-reference" },
    ],
    { xLabel: "O/F", yLabel: "Isp,vac (s)" },
  );
  renderLineChart(
    document.getElementById("of-cf-chart"),
    [
      { data: points, x: "oxidizer_fuel_ratio", y: "equilibrium_thrust_coefficient", className: "series-1" },
      { data: points, x: "oxidizer_fuel_ratio", y: "frozen_thrust_coefficient", className: "series-reference" },
    ],
    { xLabel: "O/F", yLabel: "C₍F₎" },
  );
}

function renderOffDesignNozzle(data) {
  const module = document.querySelector('[data-module="nozzle-off-design"]');
  renderResults(module, {
    ...data,
    shock_x_fraction: data.shock?.x_fraction,
    separation_status: data.separation.predicted ? "olası ayrılma" : "eşik aşılmadı",
    exit_ambient_ratio: data.separation.exit_to_ambient_pressure_ratio,
  });
  const warning = document.getElementById("off-design-warning");
  const messages = data.warnings.map((item) => item.message);
  messages.push(data.altitude_sweep.warning);
  warning.textContent = messages.join(" ");

  const pressureProfile = data.profile.map((point) => ({
    ...point,
    pressure_chamber_ratio: point.pressure_pa / data.chamber_pressure_pa,
  }));
  const maximumMach = Math.max(...data.profile.map((point) => point.mach), 1);
  const maximumPressureRatio = Math.max(
    ...pressureProfile.map((point) => point.pressure_chamber_ratio),
    1,
  );
  const shockMach = data.shock ? [
    { x_fraction: data.shock.x_fraction, mach: 0 },
    { x_fraction: data.shock.x_fraction, mach: maximumMach },
  ] : [];
  const shockPressure = data.shock ? [
    { x_fraction: data.shock.x_fraction, pressure_chamber_ratio: 0 },
    { x_fraction: data.shock.x_fraction, pressure_chamber_ratio: maximumPressureRatio },
  ] : [];
  const machSeries = [
    { data: data.profile, x: "x_fraction", y: "mach", className: "series-1" },
  ];
  const pressureSeries = [
    { data: pressureProfile, x: "x_fraction", y: "pressure_chamber_ratio", className: "series-1" },
  ];
  if (shockMach.length) {
    machSeries.push({ data: shockMach, x: "x_fraction", y: "mach", className: "series-2" });
    pressureSeries.push({ data: shockPressure, x: "x_fraction", y: "pressure_chamber_ratio", className: "series-2" });
  }
  renderLineChart(
    document.getElementById("off-design-mach-chart"),
    machSeries,
    { xLabel: "x / L", yLabel: "Mach" },
  );
  renderLineChart(
    document.getElementById("off-design-pressure-chart"),
    pressureSeries,
    { xLabel: "x / L", yLabel: "p / p₀" },
  );
  const altitude = data.altitude_sweep.points;
  renderLineChart(
    document.getElementById("altitude-thrust-chart"),
    [{ data: altitude, x: "altitude_m", y: "thrust_n", className: "series-1" }],
    { xLabel: "İrtifa (m)", yLabel: "F (N)" },
  );
  renderLineChart(
    document.getElementById("altitude-isp-chart"),
    [{ data: altitude, x: "altitude_m", y: "specific_impulse_s", className: "series-3" }],
    { xLabel: "İrtifa (m)", yLabel: "Isp (s)" },
  );
  renderLineChart(
    document.getElementById("altitude-pressure-ratio-chart"),
    [{ data: altitude, x: "altitude_m", y: "chamber_to_ambient_pressure_ratio", className: "series-reference" }],
    { xLabel: "İrtifa (m)", yLabel: "p₀ / pₐ" },
  );
}

async function loadIsentropicCharts() {
  const response = await fetch("/api/v1/isentropic-sweep", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ gamma: 1.4, minimum_mach: 0.05, maximum_mach: 5, point_count: 121 }),
  });
  const payload = await response.json();
  if (!response.ok) return;
  const points = payload.data.points;
  renderLineChart(
    document.getElementById("isentropic-ratios-chart"),
    [
      { data: points, x: "mach", y: "pressure_total_ratio", className: "series-1" },
      { data: points, x: "mach", y: "temperature_total_ratio", className: "series-2" },
      { data: points, x: "mach", y: "density_total_ratio", className: "series-3" },
    ],
    { xLabel: "Mach", yLabel: "Boyutsuz oran" },
  );
  renderLineChart(
    document.getElementById("area-mach-chart"),
    [{ data: points, x: "mach", y: "area_critical_ratio", className: "series-1" }],
    { xLabel: "Mach", yLabel: "A / A*" },
  );
}

async function calculate(form, endpoint, onSuccess = null) {
  const module = form.closest("[data-module]");
  const message = form.querySelector(".form-message");
  const button = form.querySelector("button[type='submit']");
  message.textContent = "Hesaplanıyor…";
  message.className = "form-message";
  button.disabled = true;

  try {
    const response = await fetch(endpoint, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(formPayload(form)),
    });
    const payload = await response.json();
    if (!response.ok) throw new Error(apiErrorMessage(payload, "Hesaplama tamamlanamadı."));
    renderResults(module, payload.data);
    if (onSuccess) onSuccess(payload.data);
    message.textContent = "Durum başarıyla çözüldü.";
    message.className = "form-message success";
  } catch (error) {
    message.textContent = error.message;
    message.className = "form-message error";
  } finally {
    button.disabled = false;
  }
}

const calculators = [
  ["isentropic-form", "/api/v1/isentropic"],
  ["shock-form", "/api/v1/normal-shock"],
  ["oblique-shock-form", "/api/v1/oblique-shock"],
  ["thermo-form", "/api/v1/ideal-gas"],
  ["thermal-properties-form", "/api/v1/thermal-properties", renderThermalProperties],
  ["equilibrium-sweep-form", "/api/v1/of-sweep", renderEquilibriumSweep],
  ["stagnation-form", "/api/v1/stagnation"],
  ["process-form", "/api/v1/thermo-process"],
  ["polytropic-form", "/api/v1/polytropic", renderPolytropicCharts],
  ["duct-flow-form", "/api/v1/duct-flow"],
  ["nozzle-form", "/api/v1/nozzle"],
  ["off-design-nozzle-form", "/api/v1/nozzle-off-design", renderOffDesignNozzle],
  ["geometry-form", "/api/v1/nozzle-contour", renderNozzleGeometry],
  ["loss-budget-form", "/api/v1/loss-budget", renderLossBudget],
  ["two-phase-form", "/api/v1/two-phase-band", renderTwoPhaseBand],
  ["nozzle-thermal-form", "/api/v1/nozzle-thermal", renderNozzleThermal],
  ["staging-form", "/api/v1/staging", renderStaging],
  ["engine-cycle-form", "/api/v1/engine-cycle", renderEngineCycle],
  ["uncertainty-form", "/api/v1/uncertainty", renderUncertainty],
  ["rocket-form", "/api/v1/rocket-equation"],
  ["thrust-form", "/api/v1/thrust"],
];

for (const [formId, endpoint, onSuccess] of calculators) {
  const form = document.getElementById(formId);
  form.addEventListener("submit", (event) => {
    event.preventDefault();
    calculate(form, endpoint, onSuccess);
  });
  form.requestSubmit();
}

const burnForm = document.getElementById("burn-form");
burnForm.addEventListener("submit", (event) => {
  event.preventDefault();
  runBurn(burnForm);
});
burnForm.requestSubmit();

document.getElementById("burn-json").addEventListener("click", () => {
  if (latestBurnResult) {
    downloadText(
      "rocket-burn.json",
      `${JSON.stringify(latestBurnResult, null, 2)}\n`,
      "application/json;charset=utf-8",
    );
  }
});

document.getElementById("burn-report").addEventListener("click", async () => {
  if (!latestBurnResult) return;
  const message = burnForm.querySelector(".form-message");
  try {
    const response = await fetch("/api/v1/burn/report", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ burn_result: latestBurnResult }),
    });
    const payload = await response.json();
    if (!response.ok) throw new Error(apiErrorMessage(payload, "Burn raporu üretilemedi."));
    downloadText("rocket-burn-report.html", payload.data.html, "text/html;charset=utf-8");
    message.textContent = "Burn HTML raporu indirildi.";
    message.className = "form-message success";
  } catch (error) {
    message.textContent = error.message;
    message.className = "form-message error";
  }
});

document.getElementById("burn-sidera").addEventListener("click", async () => {
  if (!latestBurnResult) return;
  const message = document.getElementById("sidera-message");
  try {
    const response = await fetch("/api/v1/burn/export/sidera", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        burn_result: latestBurnResult,
        t_start_s: document.getElementById("sidera-start").value,
        frame: document.getElementById("sidera-frame").value,
        mode: document.getElementById("sidera-mode").value,
        direction: ["sidera-x", "sidera-y", "sidera-z"].map(
          (id) => document.getElementById(id).value,
        ),
      }),
    });
    const payload = await response.json();
    if (!response.ok) throw new Error(apiErrorMessage(payload, "Sidera planı üretilemedi."));
    downloadText(
      "sidera-maneuver-plan.json",
      `${JSON.stringify(payload.data, null, 2)}\n`,
      "application/json;charset=utf-8",
    );
    message.textContent = document.getElementById("sidera-mode").value === "exact"
      ? "Exact constant-burn Sidera planı indirildi."
      : "Yaklaşık eşdeğer-sabit Sidera planı indirildi; manifest uyarılarını inceleyin.";
    message.className = "form-message success";
  } catch (error) {
    message.textContent = error.message;
    message.className = "form-message error";
  }
});

async function loadSideraCapabilities() {
  const label = document.getElementById("sidera-capability");
  try {
    const response = await fetch("/api/v1/integrations/sidera/capabilities");
    const payload = await response.json();
    if (!response.ok) throw new Error("Yetenek sorgusu başarısız.");
    const capability = payload.data;
    label.textContent = capability.constant_burn_exact
      ? `Native Sidera ${capability.version}: exact constant burn destekleniyor.`
      : "Native Sidera kullanılamıyor; doğrulanmış JSON hand-off yine üretilebilir.";
  } catch (error) {
    label.textContent = error.message;
  }
}
loadSideraCapabilities();

document.querySelectorAll("[data-geometry-export]").forEach((button) => {
  button.addEventListener("click", async () => {
    try {
      await downloadGeometryExport(button.dataset.geometryExport);
    } catch (error) {
      const message = document.querySelector("#geometry-form .form-message");
      message.textContent = error.message;
      message.className = "form-message error";
    }
  });
});

document.getElementById("workspace-save").addEventListener("click", async () => {
  const message = document.getElementById("workspace-message");
  try {
    const data = await normalizedWorkspace(loadedWorkspace || analysisWorkspace());
    loadedWorkspace = data.workspace;
    downloadText("rocket-study.rplab.json", data.normalized_json, "application/json;charset=utf-8");
    message.textContent = `${loadedWorkspace.cases.length} vaka içeren çalışma indirildi.`;
  } catch (error) {
    message.textContent = error.message;
  }
});

document.getElementById("workspace-report").addEventListener("click", async () => {
  const message = document.getElementById("workspace-message");
  try {
    const data = await normalizedWorkspace(loadedWorkspace || analysisWorkspace(), true);
    loadedWorkspace = data.workspace;
    downloadText("rocket-study-report.html", data.report.html, "text/html;charset=utf-8");
    message.textContent = "Yazdırılabilir HTML raporu indirildi.";
  } catch (error) {
    message.textContent = error.message;
  }
});

document.getElementById("workspace-open").addEventListener("change", async (event) => {
  const message = document.getElementById("workspace-message");
  const [file] = event.target.files;
  if (!file) return;
  try {
    const data = await normalizedWorkspace(JSON.parse(await file.text()));
    loadedWorkspace = data.workspace;
    const burnCase = loadedWorkspace.cases.find(
      (item) => item.model === "standalone-constant-burn" && item.results?.artifact,
    );
    if (burnCase) renderBurn(burnCase.results.artifact, burnCase.inputs);
    message.textContent = `${loadedWorkspace.name}: ${loadedWorkspace.cases.length} vaka açıldı.`;
  } catch (error) {
    message.textContent = error.message;
  }
});

document.getElementById("workspace-duplicate").addEventListener("click", async () => {
  const message = document.getElementById("workspace-message");
  try {
    const workspace = structuredClone(loadedWorkspace || analysisWorkspace());
    if (!workspace.cases.length) throw new Error("Çoğaltılacak vaka bulunamadı.");
    const copy = structuredClone(workspace.cases[0]);
    const identifiers = new Set(workspace.cases.map((item) => item.case_id));
    let suffix = 1;
    while (identifiers.has(`${copy.case_id}-copy-${suffix}`)) suffix += 1;
    copy.case_id = `${copy.case_id}-copy-${suffix}`;
    copy.name = `${copy.name} kopya ${suffix}`;
    workspace.cases.push(copy);
    const data = await normalizedWorkspace(workspace);
    loadedWorkspace = data.workspace;
    message.textContent = `${copy.name} oluşturuldu; toplam ${loadedWorkspace.cases.length} vaka.`;
  } catch (error) {
    message.textContent = error.message;
  }
});

loadIsentropicCharts();

const equilibriumPair = document.getElementById("equilibrium-pair");
equilibriumPair.addEventListener("change", () => {
  const form = document.getElementById("equilibrium-sweep-form");
  const ranges = {
    LOX_LH2: [4, 8],
    LOX_CH4: [2, 5],
  };
  const [minimum, maximum] = ranges[equilibriumPair.value];
  form.elements.minimum_ratio.value = minimum;
  form.elements.maximum_ratio.value = maximum;
});

const categoryLabels = { liquid: "Sıvı–sıvı", solid: "Katı", hybrid: "Hibrit" };

function textElement(tag, text, className) {
  const element = document.createElement(tag);
  element.textContent = text;
  if (className) element.className = className;
  return element;
}

function tradeColumn(title, items) {
  const column = document.createElement("div");
  column.append(textElement("h4", title));
  const list = document.createElement("ul");
  for (const item of items) list.append(textElement("li", item));
  column.append(list);
  return column;
}

function propellantCard(record) {
  const card = document.createElement("article");
  card.className = "propellant-card";

  const header = document.createElement("header");
  const titleBlock = document.createElement("div");
  titleBlock.append(textElement("span", record.key.toUpperCase(), "eyebrow"));
  titleBlock.append(textElement("h3", record.name));
  header.append(titleBlock, textElement("span", categoryLabels[record.category], "category-tag"));

  const isp = document.createElement("div");
  isp.className = "isp-range";
  isp.append(
    textElement("span", "TEMSİLÎ VAKUM ISP"),
    textElement("strong", `${record.vacuum_isp_min_s}–${record.vacuum_isp_max_s} s`),
  );

  const details = document.createElement("dl");
  for (const [label, value] of [
    ["Yakıt", record.fuel],
    ["Oksitleyici", record.oxidizer],
    ["Depolama", record.storage],
    ["Ateşleme", record.ignition],
    ["Karışım notu", record.mixture_ratio_note],
  ]) {
    const row = document.createElement("div");
    row.append(textElement("dt", label), textElement("dd", value));
    details.append(row);
  }

  const trades = document.createElement("div");
  trades.className = "trade-list";
  trades.append(
    tradeColumn("Güçlü yönler", record.strengths),
    tradeColumn("Sınırlamalar", record.limitations),
  );
  card.append(header, isp, details, trades);
  return card;
}

async function loadPropellants(form) {
  const grid = document.getElementById("propellant-grid");
  const message = form.querySelector(".form-message");
  const button = form.querySelector("button");
  button.disabled = true;
  message.textContent = "Katalog yükleniyor…";
  try {
    const response = await fetch("/api/v1/propellants", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(formPayload(form)),
    });
    const payload = await response.json();
    if (!response.ok) throw new Error(apiErrorMessage(payload, "Katalog yüklenemedi."));
    grid.replaceChildren(...payload.data.records.map(propellantCard));
    message.textContent = `${payload.data.count} sistem gösteriliyor.`;
    message.className = "form-message success";
  } catch (error) {
    message.textContent = error.message;
    message.className = "form-message error";
  } finally {
    button.disabled = false;
  }
}

const propellantForm = document.getElementById("propellant-form");
propellantForm.addEventListener("submit", (event) => {
  event.preventDefault();
  loadPropellants(propellantForm);
});
propellantForm.requestSubmit();

