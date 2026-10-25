/* Rocket Propulsion Lab — geçici rejim modülleri (besleme, elektrik, belirsizlik,
 * zaman çizelgesi, SASS).
 *
 * Bu dosya hesap yapmaz: her sayı /api/v1/transient/* ve /api/v1/sass/mission
 * uçlarından gelir (src/rocket_propulsion/api/transient_routes.py). Seriler
 * {x, columns} sütun tablosu olarak gelir; burada yalnız satıra çevrilip çizilir.
 * Genel işlevler app.js'ten (calculate, renderResults, renderLineChart) ve
 * visuals.js'ten (window.labCharts) gelir.                                       */

(() => {
  const byId = (id) => document.getElementById(id);

  function rows(group) {
    const names = Object.keys(group.columns);
    return Array.from({ length: group.point_count }, (_, index) => Object.fromEntries(
      names.map((name) => [name, group.columns[name][index]]),
    ));
  }

  function line(id, data, x, ys, labels) {
    renderLineChart(byId(id), ys.map(([y, slot, extra = {}]) => ({ data, x, y, className: `series-${slot}`, ...extra })), labels);
  }

  function warn(id, data) {
    const node = byId(id);
    if (node) node.textContent = (data.warnings || []).map((item) => item.message || item).join(" ");
  }

  function scoped(moduleName, values) {
    renderResults(document.querySelector(`[data-module="${moduleName}"]`), values);
  }

  function statistic(list, metric, key) {
    const item = (list || []).find((entry) => entry.metric === metric);
    return item ? item[key] : null;
  }

  const number = (value) => formatValue(value, "value");

  const renderers = {
    "transient-blowdown-form": ["/api/v1/transient/blowdown", (data) => {
      scoped("transient-blowdown", data.summary);
      const history = rows(data.series.history);
      line("blowdown-pressure-chart", history, "time_s", [["supply_pressure_pa", 1]], { xLabel: "t (s)", yLabel: "p (Pa)" });
      line("blowdown-thrust-chart", history, "time_s", [["delivered_thrust_n", 2]], { xLabel: "t (s)", yLabel: "F (N)" });
      line("blowdown-mass-chart", history, "time_s", [["remaining_propellant_mass_kg", 3]], { xLabel: "t (s)", yLabel: "m (kg)" });
      warn("blowdown-warning", data);
    }],

    "transient-feed-form": ["/api/v1/transient/regulated-feed", (data) => {
      const { propulsion, feed } = data.summary;
      scoped("transient-feed", {
        axial_thrust_n: propulsion.axial_thrust_n,
        engine_mass_flow_kg_s: propulsion.engine_mass_flow_kg_s,
        injector_pressure_pa: feed.injector_pressure_pa,
        feed_line_pressure_drop_pa: feed.feed_line_pressure_drop_pa,
        regulator_command: feed.regulator_command,
        regulator_mass_flow_kg_s: feed.regulator_mass_flow_kg_s,
        regulator_choked_label: feed.regulator_choked ? "evet" : "hayır",
      });
      line("feed-demand-chart", rows(data.series.demand_sweep), "propellant_mass_flow_kg_s",
        [["injector_pressure_pa", 1]], { xLabel: "ṁ (kg/s)", yLabel: "p (Pa)" });
      line("feed-regulator-chart", rows(data.series.regulator_sweep), "tank_pressure_pa",
        [["regulator_mass_flow_kg_s", 2]], { xLabel: "p_tank (Pa)", yLabel: "ṁ (kg/s)" });
      warn("feed-warning", data);
    }],

    "transient-line-form": ["/api/v1/transient/feed-line", (data) => {
      scoped("transient-line", data.summary.line);
      const sweep = rows(data.series.line_flow_sweep);
      line("line-damping-chart", sweep, "line_mass_flow_kg_s", [["damping_ratio", 1]], { xLabel: "ṁ (kg/s)", yLabel: "ζ" });
      line("line-friction-chart", sweep, "line_mass_flow_kg_s", [["friction_pressure_drop_pa", 2]], { xLabel: "ṁ (kg/s)", yLabel: "Δp (Pa)" });
      warn("line-warning", data);
    }],

    "transient-electric-form": ["/api/v1/transient/electric", (data) => {
      scoped("transient-electric", data.summary);
      const sweep = rows(data.series.power_sweep);
      line("electric-thrust-chart", sweep, "available_power_w", [["thrust_n", 1]], { xLabel: "P (W)", yLabel: "F (N)" });
      line("electric-level-chart", sweep, "available_power_w",
        [["level_index_increasing", 1], ["level_index_decreasing", 2]], { xLabel: "P (W)", yLabel: "Seviye" });
      line("electric-isp-chart", sweep, "available_power_w", [["specific_impulse_s", 3]], { xLabel: "P (W)", yLabel: "Isp (s)" });
      warn("electric-warning", data);
    }],

    "transient-pulse-form": ["/api/v1/transient/pulse-mode", (data) => {
      scoped("transient-pulse", {
        ...data.summary,
        pulse_count: data.summary.plan.pulse_count,
        achieved_delta_v_m_s: data.summary.plan.delivered_m_s,
      });
      const curve = rows(data.series.impulse_bit_curve);
      line("pulse-impulse-chart", curve, "on_time_s", [["impulse_n_s", 1]], { xLabel: "t_on (s)", yLabel: "I (N·s)" });
      line("pulse-isp-chart", curve, "on_time_s", [["specific_impulse_s", 3]], { xLabel: "t_on (s)", yLabel: "Isp (s)" });
      warn("pulse-warning", data);
    }],

    "transient-uncertainty-form": ["/api/v1/transient/uncertainty", (data) => {
      const stats = data.statistics;
      scoped("transient-uncertainty", {
        impulse_p50: statistic(stats, "delivered_impulse_n_s", "p50"),
        impulse_p05: statistic(stats, "delivered_impulse_n_s", "p05"),
        impulse_p95: statistic(stats, "delivered_impulse_n_s", "p95"),
        duration_p50: statistic(stats, "duration_s", "p50"),
        duration_band: `${number(statistic(stats, "duration_s", "p05"))} – ${number(statistic(stats, "duration_s", "p95"))}`,
        convergence_label: data.convergence.converged
          ? `yakınsadı (ortalama değişimi ${number(data.convergence.mean_relative_change)})`
          : "yakınsamadı — örnek sayısını artır",
      });
      const realizations = rows(data.series.realizations);
      const charts = window.labCharts;
      charts.renderHistogram(byId("uncertainty-impulse-chart"), realizations.map((row) => row.delivered_impulse_n_s), {
        xLabel: "Toplam impuls (N·s)",
        references: [{ x: statistic(stats, "delivered_impulse_n_s", "p50"), label: "P50" }],
      });
      charts.renderHistogram(byId("uncertainty-duration-chart"), realizations.map((row) => row.duration_s), {
        xLabel: "Yanış süresi (s)",
        references: [{ x: statistic(stats, "duration_s", "p50"), label: "P50" }],
      });
      charts.renderBars(byId("uncertainty-sensitivity-chart"), data.sensitivity
        .filter((item) => item.metric === "delivered_impulse_n_s" && item.correlation !== null)
        .map((item) => ({ label: item.parameter, segments: [{ value: item.correlation, slot: 1, name: "korelasyon" }] })),
      { unit: "", xLabel: "Pearson korelasyonu" });
      warn("uncertainty-warning", data);
    }],

    "transient-scenarios-form": ["/api/v1/transient/scenarios", (data) => {
      scoped("transient-scenarios", {
        ...data.summary,
        expected_impulse: statistic(data.statistics, "delivered_impulse_n_s", "mean"),
      });
      const table = rows(data.series.consequences);
      const charts = window.labCharts;
      charts.renderBars(byId("scenario-impulse-chart"), table.map((row) => ({
        label: row.scenario_id, segments: [{ value: row.delivered_impulse_n_s, slot: 1, name: "impuls" }],
      })), { unit: "N·s", xLabel: "Toplam impuls (N·s)" });
      charts.renderBars(byId("scenario-probability-chart"), table.map((row) => ({
        label: row.scenario_id, segments: [{ value: row.probability, slot: 4, name: "olasılık" }],
      })), { unit: "", xLabel: "Olasılık" });
      warn("scenario-warning", data);
    }],

    "transient-calibration-form": ["/api/v1/transient/calibration", (data) => {
      scoped("transient-calibration", {
        ...data.summary,
        observation_label: `${data.summary.observation_count} / ${data.summary.engine_count}`,
      });
      window.labCharts.renderRanges(byId("calibration-engines-chart"), rows(data.series.engines).map((row) => ({
        label: row.engine_id, min: row.minimum_scale, max: row.maximum_scale, slot: 1, group: "koşu aralığı",
      })), { xLabel: "Ölçek", unit: "" });
      warn("calibration-warning", data);
    }],

    "transient-timeline-form": ["/api/v1/transient/timeline", (data) => {
      scoped("transient-timeline", {
        ...data.summary,
        event_label: `${data.summary.burn_count} / ${data.summary.jettison_count}`,
      });
      const history = rows(data.series.history);
      line("timeline-thrust-chart", history, "time_s", [["axial_thrust_n", 1]], { xLabel: "t (s)", yLabel: "F (N)" });
      line("timeline-flow-chart", history, "time_s", [["mass_flow_kg_s", 2]], { xLabel: "t (s)", yLabel: "ṁ (kg/s)" });
      line("timeline-jettison-chart", history, "time_s", [["jettisoned_mass_kg", 4]], { xLabel: "t (s)", yLabel: "m (kg)" });
      warn("timeline-warning", data);
    }],

    "transient-sass-form": ["/api/v1/sass/mission", (data) => {
      const scored = data.candidates.filter((item) => item.score !== null);
      const best = scored.reduce((top, item) => (top === null || item.score > top.score ? item : top), null);
      scoped("transient-sass", {
        ...data.summary,
        best_candidate: best ? best.identifier : "uygun seçenek yok",
        eligible_label: `${data.summary.eligible_count} / ${data.summary.candidate_count}`,
      });
      const trade = rows(data.series.trade);
      const charts = window.labCharts;
      charts.renderBars(byId("sass-mass-chart"), trade.filter((row) => row.wet_mass_kg !== null).map((row) => ({
        label: row.identifier, segments: [{ value: row.wet_mass_kg, slot: 2, name: "ıslak kütle" }],
      })), { unit: "kg", xLabel: "Islak kütle (kg)" });
      charts.renderBars(byId("sass-score-chart"), trade.filter((row) => row.score !== null).map((row) => ({
        label: row.identifier, segments: [{ value: row.score, slot: 1, name: "puan" }],
      })), { unit: "", xLabel: "Ağırlıklı puan (0–1)" });
      warn("sass-warning", data);
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
