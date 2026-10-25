/* Rocket Propulsion Lab — grafik etkileşimi ve ek görseller.
 *
 * Bu katman hesap yapmaz: her sayı Python API'sinden gelir. Taramalar mevcut uçları
 * (normal-shock, oblique-shock, duct-flow, rocket-equation, nozzle-contour) farklı
 * girdilerle çağırır; tarayıcıda yalnız ölçekleme, sayma (histogram) ve çizim yapılır.
 *
 * Sunumatik desenleri (vendor/sunumatik, PROVENANCE.json):
 *   · chart-theme.css seri yuvaları (--color-data-N) ve kesik çizgi ikincil kodlaması
 *   · hover_effects "m-ipucu" (değer ipucu), "k-odak" (odak dışını soldurma),
 *     "k-kaldir" (kart kaldırma) — ölçüleri bu uygulamaya uyarlandı, ≤ 300 ms
 *   · grafik-cekirdek 1–2–2,5–5 tik kuralı (app.js axisTicks) ve tipografik eksi
 * Genel işlevler app.js'ten gelir: formatValue, svgElement, axisTicks, chartTick.     */

(() => {
  const W = 720;
  const H = 280;

  /* ── yardımcılar ───────────────────────────────────────────────────── */

  const el = (tag, attributes = {}, text = "") => svgElement(tag, attributes, text);
  const fmt = (value) => formatValue(value, "value");
  // On-chart labels: four significant digits (datawiz: no excess digits); tooltips keep fmt.
  const shortNumber = new Intl.NumberFormat("tr-TR", { maximumSignificantDigits: 4 });
  const short = (value) => shortNumber.format(value);
  const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));

  async function post(url, body) {
    const response = await fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    const payload = await response.json();
    if (!response.ok) throw new Error(apiErrorMessage(payload, "İstek başarısız."));
    return payload.data;
  }

  /** Run promise factories with bounded concurrency; failed items resolve to null. */
  async function pool(factories, limit = 6) {
    const results = new Array(factories.length).fill(null);
    let next = 0;
    async function worker() {
      while (next < factories.length) {
        const index = next++;
        try { results[index] = await factories[index](); } catch { results[index] = null; }
      }
    }
    await Promise.all(Array.from({ length: Math.min(limit, factories.length) }, worker));
    return results;
  }

  const cache = new Map();
  function cached(key, factory) {
    if (!cache.has(key)) cache.set(key, factory().catch((error) => { cache.delete(key); throw error; }));
    return cache.get(key);
  }

  const byId = (id) => document.getElementById(id);
  const values = (form) => Object.fromEntries(new FormData(form).entries());

  /** Linear axis object compatible with app.js axisTicks (fraction ↔ value). */
  function linearAxis(lower, upper) {
    const span = upper - lower || 1;
    return { tick: (fraction) => lower + fraction * span, fraction: (v) => (v - lower) / span };
  }

  /* ── ipucu (m-ipucu) ───────────────────────────────────────────────── */

  const tip = document.createElement("div");
  tip.className = "chart-tip";
  tip.setAttribute("role", "status");
  tip.hidden = true;
  document.body.append(tip);

  function showTip(html, clientX, clientY) {
    tip.innerHTML = html;
    tip.hidden = false;
    const box = tip.getBoundingClientRect();
    const x = Math.min(window.innerWidth - box.width - 12, clientX + 14);
    const y = clientY - box.height - 14 < 8 ? clientY + 18 : clientY - box.height - 14;
    tip.style.transform = `translate(${Math.max(8, x)}px, ${y}px)`;
  }
  const hideTip = () => { tip.hidden = true; };
  const escape = (text) => String(text).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);

  function viewBoxPoint(svg, event) {
    const matrix = svg.getScreenCTM();
    if (!matrix) return null;
    const point = new DOMPoint(event.clientX, event.clientY).matrixTransform(matrix.inverse());
    return { x: point.x, y: point.y };
  }

  /* ── çizgi grafiklerinde imleç okuması ve seri vurgusu ─────────────── */

  function legendNames(svg) {
    return [...(svg.closest("figure")?.querySelectorAll(".chart-legend span") || [])].map((s) => s.textContent);
  }

  function bindLineHover(svg) {
    if (svg.__hoverBound) return;
    svg.__hoverBound = true;
    svg.addEventListener("pointermove", (event) => {
      const chart = svg.__chart;
      if (!chart) return;
      const p = viewBoxPoint(svg, event);
      const { left, top, width, height } = chart.plot;
      svg.querySelector(".chart-hover")?.remove();
      if (!p || p.x < left || p.x > left + width || p.y < top || p.y > top + height) { hideTip(); return; }
      const names = legendNames(svg);
      const group = el("g", { class: "chart-hover chart-deco" });
      const rows = [];
      let anchorX = null;
      chart.series.forEach((series, index) => {
        let best = null;
        for (const point of series.points) {
          const d = Math.abs(chart.xScale(point.x) - p.x);
          if (!best || d < best.d) best = { d, point };
        }
        if (!best || best.d > width / 8) return;
        const cx = chart.xScale(best.point.x);
        const cy = chart.yScale(best.point.y);
        anchorX ??= { cx, x: best.point.x };
        group.append(el("circle", { cx: cx.toFixed(1), cy: cy.toFixed(1), r: 4.5, class: `hover-dot ${series.className}`, "data-series-slot": series.className.replace("series-", "") }));
        rows.push(`<span class="tip-row"><i class="tip-swatch ${series.className}"></i>${escape(names[index] || chart.yLabel)}<b>${fmt(best.point.y)}</b></span>`);
      });
      if (!anchorX) { hideTip(); return; }
      group.prepend(el("line", { x1: anchorX.cx, x2: anchorX.cx, y1: top, y2: top + height, class: "hover-guide" }));
      svg.append(group);
      showTip(`<span class="tip-head">${escape(chart.xLabel)} = ${fmt(anchorX.x)}</span>${rows.join("")}`, event.clientX, event.clientY);
    });
    svg.addEventListener("pointerleave", () => { svg.querySelector(".chart-hover")?.remove(); hideTip(); });
  }

  // k-odak: lejanttaki bir seri üzerine gelince diğerleri solar.
  document.addEventListener("pointerover", (event) => {
    const item = event.target.closest?.(".chart-legend span");
    if (!item) return;
    const figure = item.closest("figure");
    const index = [...item.parentElement.children].indexOf(item);
    figure.querySelectorAll(".series-line").forEach((line, i) => line.classList.toggle("is-dimmed", i !== index));
  });
  document.addEventListener("pointerout", (event) => {
    const item = event.target.closest?.(".chart-legend span");
    if (item) item.closest("figure").querySelectorAll(".is-dimmed").forEach((line) => line.classList.remove("is-dimmed"));
  });

  // Çubuk, nokta ve alan şeritleri: data-tip taşıyan her öğe ipucu verir.
  document.addEventListener("pointermove", (event) => {
    const target = event.target.closest?.("[data-tip]");
    if (target && target.closest(".sci-chart")) showTip(target.dataset.tip, event.clientX, event.clientY);
  });
  document.addEventListener("pointerout", (event) => {
    if (event.target.closest?.("[data-tip]")) hideTip();
  });

  /* ── çözülen durum işaretçisi ──────────────────────────────────────── */

  function setMarker(id, marker) {
    const svg = byId(id);
    if (!svg) return;
    svg.__marker = marker;
    drawMarker(svg);
  }

  function drawMarker(svg) {
    svg.querySelector(".chart-marker")?.remove();
    const chart = svg.__chart;
    const marker = svg.__marker;
    if (!chart || !marker || !Number.isFinite(marker.x)) return;
    const { top, height, left, width } = chart.plot;
    const cx = chart.xScale(marker.x);
    if (cx < left - 1 || cx > left + width + 1) return;
    const group = el("g", { class: "chart-marker chart-deco" });
    group.append(el("line", { x1: cx, x2: cx, y1: top, y2: top + height, class: "marker-line" }));
    for (const y of marker.ys || []) {
      if (Number.isFinite(y)) group.append(el("circle", { cx, cy: chart.yScale(y), r: 5.5, class: "marker-dot" }));
    }
    const anchor = cx > left + width * 0.75 ? "end" : "start";
    group.append(el("text", { x: cx + (anchor === "end" ? -6 : 6), y: top + 13, "text-anchor": anchor, class: "marker-label" }, marker.label));
    svg.append(group);
  }

  /* ── yeni çizimlerde efektler: beliriş + etkileşim + işaretçi ──────── */

  function decorate(svg) {
    if (svg.__chart) { bindLineHover(svg); drawMarker(svg); }
    if (!svg.__revealed) {
      svg.__revealed = true;
      svg.classList.add("is-revealing");
      setTimeout(() => svg.classList.remove("is-revealing"), 320);
    }
  }

  new MutationObserver((records) => {
    const touched = new Set();
    for (const record of records) {
      const svg = record.target;
      if (!(svg instanceof SVGSVGElement) || !svg.classList.contains("sci-chart")) continue;
      // Only a real redraw counts. Adding or removing our own hover/marker layers must not
      // re-trigger decorate(), or marker redraws would feed this observer forever.
      const nodes = [...record.addedNodes, ...record.removedNodes];
      if (nodes.every((node) => node.classList?.contains("chart-deco"))) continue;
      touched.add(svg);
    }
    touched.forEach(decorate);
  }).observe(document.querySelector("main"), { childList: true, subtree: true });

  /* ── genel SVG çizim parçaları ─────────────────────────────────────── */

  function frame(svg) {
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    svg.replaceChildren();
  }

  function gridAndAxes(svg, plot, xAxis, yAxis, { xLabel, yLabel, xTicks = true, yTicks = true, integerX = false }) {
    if (xTicks) {
      for (const { fraction, label } of axisTicks(xAxis, { integer: integerX })) {
        const x = plot.left + fraction * plot.width;
        svg.append(el("line", { x1: x, x2: x, y1: plot.top, y2: plot.top + plot.height, class: "chart-grid-line" }),
          el("text", { x, y: plot.top + plot.height + 20, "text-anchor": "middle", class: "chart-tick" }, label));
      }
    }
    if (yTicks) {
      for (const { fraction, label } of axisTicks(yAxis)) {
        const y = plot.top + (1 - fraction) * plot.height;
        svg.append(el("line", { x1: plot.left, x2: plot.left + plot.width, y1: y, y2: y, class: "chart-grid-line" }),
          el("text", { x: plot.left - 10, y: y + 4, "text-anchor": "end", class: "chart-tick" }, label));
      }
    }
    svg.append(
      el("line", { x1: plot.left, x2: plot.left + plot.width, y1: plot.top + plot.height, y2: plot.top + plot.height, class: "chart-axis" }),
      el("line", { x1: plot.left, x2: plot.left, y1: plot.top, y2: plot.top + plot.height, class: "chart-axis" }),
    );
    if (xLabel) svg.append(el("text", { x: plot.left + plot.width / 2, y: H - 8, "text-anchor": "middle", class: "chart-axis-title" }, xLabel));
    if (yLabel) svg.append(el("text", { x: 16, y: plot.top + plot.height / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${plot.top + plot.height / 2})`, class: "chart-axis-title" }, yLabel));
  }

  /**
   * Horizontal bars from zero. rows: [{label, segments: [{value, slot, name}]}].
   * Stacked segments share a row; negative values extend left of zero.
   */
  function renderBars(svg, rows, { unit, xLabel }) {
    frame(svg);
    if (!rows.length) return;
    // Stacked segments carry meaning by colour; name them on the chart so colour is not the
    // only cue (datawiz: renk tek ipucu olamaz).
    const kinds = new Map();
    for (const row of rows) for (const s of row.segments) if (s.name && !kinds.has(s.slot)) kinds.set(s.slot, s.name);
    const stacked = rows.some((row) => row.segments.length > 1);
    const top = stacked ? 34 : 12;
    const plot = { left: 168, top, width: W - 168 - 64, height: H - top - 48 };
    if (stacked) {
      let x = plot.left;
      for (const [slot, name] of kinds) {
        svg.append(el("rect", { x, y: 8, width: 14, height: 14, rx: 3, class: `bar series-${slot} legend-key`, "data-series-slot": slot }),
          el("text", { x: x + 20, y: 19, class: "bar-label" }, name));
        x += 34 + name.length * 7.5;
      }
    }
    const totals = rows.map((row) => row.segments.reduce((sum, s) => sum + Math.max(0, s.value), 0));
    const negatives = rows.map((row) => row.segments.reduce((sum, s) => sum + Math.min(0, s.value), 0));
    const xMax = Math.max(...totals, 0);
    const xMin = Math.min(...negatives, 0);
    // headroom for the end-of-bar value label
    const xAxis = linearAxis(xMin, xMax === xMin ? xMin + 1 : xMax * 1.2);
    gridAndAxes(svg, plot, xAxis, null, { xLabel, yTicks: false });
    const band = plot.height / rows.length;
    const barHeight = Math.min(34, band * 0.62);
    const zero = plot.left + xAxis.fraction(0) * plot.width;
    rows.forEach((row, index) => {
      const cy = plot.top + band * (index + 0.5);
      svg.append(el("text", { x: plot.left - 12, y: cy + 4, "text-anchor": "end", class: "bar-label" }, row.label));
      let positive = 0;
      let negative = 0;
      for (const segment of row.segments) {
        const start = segment.value >= 0 ? positive : negative + segment.value;
        if (segment.value >= 0) positive += segment.value; else negative += segment.value;
        const x0 = plot.left + xAxis.fraction(start) * plot.width;
        const width = Math.abs(segment.value) / (xAxis.tick(1) - xAxis.tick(0)) * plot.width;
        svg.append(el("rect", {
          x: x0.toFixed(1), y: (cy - barHeight / 2).toFixed(1), width: Math.max(0.5, width).toFixed(1), height: barHeight.toFixed(1),
          rx: 3, class: `bar series-${segment.slot}`, "data-series-slot": segment.slot,
          "data-tip": `<span class="tip-head">${escape(row.label)}</span><span class="tip-row"><i class="tip-swatch series-${segment.slot}"></i>${escape(segment.name || "")}<b>${fmt(segment.value)} ${escape(unit)}</b></span>`,
        }));
      }
      const end = plot.left + xAxis.fraction(positive) * plot.width;
      svg.append(el("text", { x: end + 6, y: cy + 4, class: "bar-value" }, `${short(positive + negative)} ${unit}`.trim()));
    });
    svg.append(el("line", { x1: zero, x2: zero, y1: plot.top, y2: plot.top + plot.height, class: "chart-axis" }));
  }

  /** Vertical histogram from zero with optional reference lines [{x, label}]. */
  function renderHistogram(svg, samples, { xLabel, references = [] }) {
    frame(svg);
    const finite = samples.filter(Number.isFinite);
    if (finite.length < 2) return;
    const min = Math.min(...finite);
    const max = Math.max(...finite);
    const bins = Math.max(5, Math.min(30, Math.ceil(Math.sqrt(finite.length))));
    const width = (max - min) / bins || 1;
    const counts = new Array(bins).fill(0);
    for (const value of finite) counts[Math.min(bins - 1, Math.floor((value - min) / width))] += 1;
    const plot = { left: 72, top: 16, width: W - 72 - 20, height: H - 16 - 48 };
    const xAxis = linearAxis(min, max);
    const yAxis = linearAxis(0, Math.max(...counts) * 1.08);
    gridAndAxes(svg, plot, xAxis, yAxis, { xLabel, yLabel: "Örnek sayısı" });
    counts.forEach((count, index) => {
      const x0 = plot.left + index / bins * plot.width;
      const h = count / yAxis.tick(1) * plot.height;
      svg.append(el("rect", {
        x: (x0 + 1).toFixed(1), y: (plot.top + plot.height - h).toFixed(1), width: Math.max(1, plot.width / bins - 2).toFixed(1), height: h.toFixed(1),
        class: "bar series-1", "data-series-slot": "1",
        "data-tip": `<span class="tip-head">${fmt(min + index * width)} – ${fmt(min + (index + 1) * width)}</span><span class="tip-row">Örnek<b>${count}</b></span>`,
      }));
    });
    for (const reference of references) {
      const x = plot.left + xAxis.fraction(reference.x) * plot.width;
      svg.append(el("line", { x1: x, x2: x, y1: plot.top, y2: plot.top + plot.height, class: "marker-line" }),
        el("text", { x: x + 5, y: plot.top + 12, class: "marker-label" }, reference.label));
    }
  }

  /** Min–max ranges per row (dumbbell), slot per category. */
  function renderRanges(svg, rows, { xLabel, unit }) {
    frame(svg);
    if (!rows.length) return;
    const plot = { left: 168, top: 12, width: W - 168 - 70, height: H - 12 - 48 };
    const low = Math.min(...rows.map((r) => r.min));
    const high = Math.max(...rows.map((r) => r.max));
    const pad = (high - low) * 0.06 || 1;
    const xAxis = linearAxis(low - pad, high + pad);
    gridAndAxes(svg, plot, xAxis, null, { xLabel, yTicks: false });
    const band = plot.height / rows.length;
    rows.forEach((row, index) => {
      const cy = plot.top + band * (index + 0.5);
      const x0 = plot.left + xAxis.fraction(row.min) * plot.width;
      const x1 = plot.left + xAxis.fraction(row.max) * plot.width;
      const tipText = `<span class="tip-head">${escape(row.label)}</span><span class="tip-row"><i class="tip-swatch series-${row.slot}"></i>${escape(row.group)}<b>${fmt(row.min)}–${fmt(row.max)} ${escape(unit)}</b></span>`;
      svg.append(
        el("text", { x: plot.left - 12, y: cy + 4, "text-anchor": "end", class: "bar-label" }, row.label),
        el("line", { x1: x0, x2: x1, y1: cy, y2: cy, class: `range-line series-${row.slot}`, "data-series-slot": row.slot, "data-tip": tipText }),
        el("circle", { cx: x0, cy, r: 6, class: `range-end series-${row.slot}`, "data-series-slot": row.slot, "data-tip": tipText }),
        el("circle", { cx: x1, cy, r: 6, class: `range-end series-${row.slot}`, "data-series-slot": row.slot, "data-tip": tipText }),
        el("text", { x: x1 + 10, y: cy + 4, class: "bar-value" }, `${short(row.min)}–${short(row.max)}`),
      );
    });
  }

  /**
   * Axisymmetric flow field: the contour drawn at true aspect ratio, each station
   * interval filled by its mean Mach number on the palette's sequential ramp.
   */
  function renderMachField(svg, stations, { caption } = {}) {
    frame(svg);
    const valid = stations.filter((s) => Number.isFinite(s.x_m) && Number.isFinite(s.radius_m) && Number.isFinite(s.mach));
    if (valid.length < 2) return;
    const xMin = Math.min(...valid.map((s) => s.x_m));
    const xMax = Math.max(...valid.map((s) => s.x_m));
    const rMax = Math.max(...valid.map((s) => s.radius_m));
    const machMin = Math.min(...valid.map((s) => s.mach));
    const machMax = Math.max(...valid.map((s) => s.mach));
    const margin = { left: 24, right: 150, top: 18, bottom: 18 };
    const scale = Math.min((W - margin.left - margin.right) / (xMax - xMin), (H - margin.top - margin.bottom) / (2 * rMax));
    const cx = margin.left + ((W - margin.left - margin.right) - (xMax - xMin) * scale) / 2;
    const cy = H / 2;
    const X = (x) => cx + (x - xMin) * scale;
    const R = (r) => r * scale;
    const share = (mach) => Math.round(100 * (mach - machMin) / (machMax - machMin || 1));
    for (let i = 0; i < valid.length - 1; i += 1) {
      const a = valid[i];
      const b = valid[i + 1];
      const mach = (a.mach + b.mach) / 2;
      svg.append(el("polygon", {
        points: [[X(a.x_m), cy - R(a.radius_m)], [X(b.x_m), cy - R(b.radius_m)], [X(b.x_m), cy + R(b.radius_m)], [X(a.x_m), cy + R(a.radius_m)]]
          .map(([x, y]) => `${x.toFixed(1)},${y.toFixed(1)}`).join(" "),
        class: "field-strip",
        // aynı renkte ince kontur: komşu şeritler arasında kenar yumuşatma çizgisi kalmasın
        style: `--strip: color-mix(in oklab, var(--ramp-seq-5, var(--color-accent)) ${share(mach)}%, var(--ramp-seq-2, var(--color-rule))); fill: var(--strip); stroke: var(--strip); stroke-width: 0.8`,
        "data-tip": `<span class="tip-head">x = ${fmt(a.x_m)} m</span><span class="tip-row">Mach<b>${fmt(mach)}</b></span><span class="tip-row">r<b>${fmt(a.radius_m)} m</b></span><span class="tip-row">p/p₀<b>${fmt(a.pressure_total_ratio)}</b></span>`,
      }));
    }
    const wall = (sign) => valid.map((s, i) => `${i ? "L" : "M"}${X(s.x_m).toFixed(1)},${(cy + sign * R(s.radius_m)).toFixed(1)}`).join(" ");
    svg.append(
      el("path", { d: wall(-1), class: "field-wall" }),
      el("path", { d: wall(1), class: "field-wall" }),
      el("line", { x1: X(xMin), x2: X(xMax), y1: cy, y2: cy, class: "field-axis" }),
    );
    const throat = valid.reduce((best, s) => (s.radius_m < best.radius_m ? s : best));
    svg.append(el("text", { x: X(throat.x_m), y: cy - R(throat.radius_m) - 8, "text-anchor": "middle", class: "field-label" }, "boğaz · M = 1"));
    const exit = valid.at(-1);
    svg.append(el("text", { x: X(exit.x_m) + 8, y: cy + 4, class: "field-label" }, `Mₑ = ${fmt(exit.mach)}`));
    // renk ölçeği: 6 basamak, doğrudan etiketli
    const legendX = W - 118;
    svg.append(el("text", { x: legendX, y: 30, class: "field-label" }, caption || "Yerel Mach"));
    for (let i = 0; i < 6; i += 1) {
      const mach = machMin + (machMax - machMin) * (5 - i) / 5;
      svg.append(
        el("rect", { x: legendX, y: 42 + i * 26, width: 22, height: 22, rx: 3, class: "field-strip", style: `fill: color-mix(in oklab, var(--ramp-seq-5, var(--color-accent)) ${share(mach)}%, var(--ramp-seq-2, var(--color-rule)))` }),
        el("text", { x: legendX + 30, y: 57 + i * 26, class: "chart-tick" }, fmt(Number(mach.toPrecision(3)))),
      );
    }
  }

  /** Line chart via app.js renderLineChart plus a state marker. */
  function lineChart(id, series, labels, marker) {
    const svg = byId(id);
    if (!svg) return;
    svg.__marker = marker || null;
    renderLineChart(svg, series, labels);
  }

  /* ── modüllere özel görseller ──────────────────────────────────────── */

  const handlers = {
    "isentropic-form"(data) {
      const label = `M = ${fmt(data.mach)}`;
      setMarker("isentropic-ratios-chart", { x: data.mach, ys: [data.pressure_total_ratio, data.temperature_total_ratio, data.density_total_ratio], label });
      setMarker("area-mach-chart", { x: data.mach, ys: [data.area_critical_ratio], label });
    },

    async "shock-form"(data, form) {
      const gamma = Number(values(form).gamma);
      const sweep = await cached(`normal:${gamma}`, async () => {
        const machs = linspace(1.02, 5, 40);
        const rows = await pool(machs.map((m) => () => post("/api/v1/normal-shock", { upstream_mach: m, gamma })));
        return rows.map((row, i) => row && { m1: machs[i], ...row }).filter(Boolean);
      });
      const m1 = Number(values(form).upstream_mach);
      lineChart("normal-shock-jump-chart", [
        { data: sweep, x: "m1", y: "pressure_ratio", className: "series-1" },
        { data: sweep, x: "m1", y: "temperature_ratio", className: "series-2" },
        { data: sweep, x: "m1", y: "density_ratio", className: "series-3" },
      ], { xLabel: "M₁", yLabel: "Sıçrama oranı" }, { x: m1, ys: [data.pressure_ratio, data.temperature_ratio, data.density_ratio], label: `M₁ = ${fmt(m1)}` });
      lineChart("normal-shock-loss-chart", [
        { data: sweep, x: "m1", y: "downstream_mach", className: "series-1" },
        { data: sweep, x: "m1", y: "stagnation_pressure_ratio", className: "series-2" },
      ], { xLabel: "M₁", yLabel: "Boyutsuz" }, { x: m1, ys: [data.downstream_mach, data.stagnation_pressure_ratio], label: `M₁ = ${fmt(m1)}` });
    },

    async "oblique-shock-form"(data, form) {
      const input = values(form);
      const m1 = Number(input.upstream_mach);
      const gamma = Number(input.gamma);
      const thetaMax = data.maximum_turning_angle_deg;
      const sweep = await cached(`oblique:${m1}:${gamma}`, async () => {
        const thetas = linspace(0.2, thetaMax * 0.999, 32);
        const call = (branch) => pool(thetas.map((t) => () => post("/api/v1/oblique-shock", { upstream_mach: m1, turning_angle_deg: t, branch, gamma })));
        const [weak, strong] = await Promise.all([call("weak"), call("strong")]);
        return {
          weak: weak.map((row, i) => row && { theta: thetas[i], beta: row.wave_angle_deg }).filter(Boolean),
          strong: strong.map((row, i) => row && { theta: thetas[i], beta: row.wave_angle_deg }).filter(Boolean),
        };
      });
      const theta = Number(input.turning_angle_deg);
      lineChart("oblique-shock-polar-chart", [
        { data: sweep.weak, x: "theta", y: "beta", className: "series-1" },
        { data: sweep.strong, x: "theta", y: "beta", className: "series-2" },
      ], { xLabel: "Akış sapması θ (°)", yLabel: "Şok açısı β (°)" }, { x: theta, ys: [data.wave_angle_deg], label: `θ = ${fmt(theta)}°` });
    },

    async "duct-flow-form"(data, form) {
      const input = values(form);
      const gamma = Number(input.gamma);
      const sweep = await cached(`duct:${input.model}:${gamma}`, async () => {
        const machs = [...linspace(0.12, 0.98, 22), ...linspace(1.02, 3, 18)];
        const rows = await pool(machs.map((m) => () => post("/api/v1/duct-flow", { model: input.model, mach: m, gamma })));
        return rows.map((row, i) => row && { mach: machs[i], ...row }).filter(Boolean);
      });
      const mach = Number(input.mach);
      lineChart("duct-flow-chart", [
        { data: sweep, x: "mach", y: "temperature_critical_ratio", className: "series-1" },
        { data: sweep, x: "mach", y: "pressure_critical_ratio", className: "series-2" },
        { data: sweep, x: "mach", y: "stagnation_pressure_critical_ratio", className: "series-3" },
      ], { xLabel: "Mach", yLabel: "Sonik referansa oran" }, { x: mach, ys: [data.temperature_critical_ratio, data.pressure_critical_ratio, data.stagnation_pressure_critical_ratio], label: `M = ${fmt(mach)}` });
      const caption = byId("duct-flow-chart")?.closest("figure")?.querySelector("figcaption");
      if (caption) caption.textContent = input.model === "fanno" ? "Fanno çizgisi" : "Rayleigh çizgisi";
    },

    async "rocket-form"(data) {
      const ve = data.effective_exhaust_velocity_m_s;
      const sweep = await cached(`rocket:${ve}`, async () => {
        const ratios = linspace(1.05, 12, 36);
        const rows = await pool(ratios.map((ratio) => () => post("/api/v1/rocket-equation", { initial_mass_kg: ratio, final_mass_kg: 1, effective_exhaust_velocity_m_s: ve })));
        return rows.map((row) => row && { ratio: row.mass_ratio, dv: row.delta_v_m_s }).filter(Boolean);
      });
      lineChart("rocket-dv-chart", [{ data: sweep, x: "ratio", y: "dv", className: "series-1" }],
        { xLabel: "Kütle oranı m₀/mₓ", yLabel: "İdeal Δv (m/s)" }, { x: data.mass_ratio, ys: [data.delta_v_m_s], label: `m₀/mₓ = ${fmt(data.mass_ratio)}` });
    },

    "thrust-form"(data) {
      renderBars(byId("thrust-parts-chart"), [
        { label: "Momentum itkisi", segments: [{ value: data.momentum_thrust_n, slot: 1, name: "ṁ·Vₑ" }] },
        { label: "Basınç itkisi", segments: [{ value: data.pressure_thrust_n, slot: 2, name: "(pₑ − pₐ)·Aₑ" }] },
        { label: "Toplam itki", segments: [{ value: data.total_thrust_n, slot: 3, name: "F" }] },
      ], { unit: "N", xLabel: "İtki (N)" });
    },

    async "nozzle-form"(data, form) {
      const input = values(form);
      const geometry = await cached(`field:${input.area_ratio}:${input.gamma}:${input.throat_area_m2}`, () => post("/api/v1/nozzle-contour", {
        contour: "rao", throat_area_m2: input.throat_area_m2, area_ratio: input.area_ratio, gamma: input.gamma,
        half_angle_deg: 15, chamber_radius_ratio: 2.5, bell_length_fraction: 0.8, characteristic_count: 20,
        throat_angle_deg: 30, exit_angle_deg: 8, station_count: 121, axisymmetric_correction: true,
        include_exports: false, include_comparison: false, angular_segments: 32,
      }));
      renderMachField(byId("nozzle-field-chart"), geometry.stations);
    },

    "geometry-form"(data) {
      renderMachField(byId("geometry-field-chart"), data.stations);
    },

    "equilibrium-sweep-form"(data) {
      const peak = data.points.reduce((best, p) => (p.equilibrium_vacuum_specific_impulse_s > best.equilibrium_vacuum_specific_impulse_s ? p : best));
      const label = `en yüksek Isp · O/F = ${fmt(peak.oxidizer_fuel_ratio)}`;
      setMarker("of-isp-chart", { x: peak.oxidizer_fuel_ratio, ys: [peak.equilibrium_vacuum_specific_impulse_s], label });
      setMarker("of-temperature-chart", { x: peak.oxidizer_fuel_ratio, ys: [peak.adiabatic_flame_temperature_k], label });
    },

    "loss-budget-form"(data) {
      const terms = data.terms.filter((t) => t.enabled).sort((a, b) => b.thrust_loss_n - a.thrust_loss_n);
      renderBars(byId("loss-terms-chart"), terms.map((t) => ({ label: t.label, segments: [{ value: t.thrust_loss_n, slot: 2, name: "İtki kaybı" }] })),
        { unit: "N", xLabel: "Terime göre itki kaybı (N)" });
    },

    "staging-form"(data) {
      const rows = data.stages.map((stage, i) => ({
        label: `${i + 1}. kademe · ${stage.name}`,
        segments: [
          { value: stage.propellant_mass_kg, slot: 1, name: "İtici" },
          { value: stage.dry_mass_kg, slot: 2, name: "Kuru yapı" },
        ],
      }));
      rows.push({ label: "Faydalı yük", segments: [{ value: data.payload_mass_kg, slot: 3, name: "Faydalı yük" }] });
      renderBars(byId("staging-mass-chart"), rows, { unit: "kg", xLabel: "Kütle (kg)" });
    },

    "engine-cycle-form"(data) {
      renderBars(byId("cycle-power-chart"), [
        { label: "Pompa gücü", segments: [{ value: data.pump_power_w, slot: 1, name: "Gerekli" }] },
        { label: "Türbin gücü", segments: [{ value: data.turbine_power_w, slot: 2, name: "Sağlanan" }] },
      ], { unit: "W", xLabel: "Mil gücü (W)" });
    },

    "uncertainty-form"(data) {
      renderHistogram(byId("uncertainty-histogram-chart"), data.samples.map((s) => s.output), {
        xLabel: "Δv (m/s)",
        references: [
          { x: data.output_percentiles.p05, label: "P05" },
          { x: data.output_percentiles.p95, label: "P95" },
        ],
      });
    },

    "burn-form"(data) {
      const final = Object.fromEntries(data.final_state.tank_masses_kg);
      const names = { oxidizer: "Oksitleyici", fuel: "Yakıt" };
      renderBars(byId("burn-tank-chart"), data.summary.consumed_by_tank_kg.map(([tank, consumed]) => ({
        label: names[tank] || tank,
        segments: [
          { value: consumed, slot: 2, name: "Tüketilen" },
          { value: final[tank], slot: 1, name: "Kalan" },
        ],
      })), { unit: "kg", xLabel: "Tank kütlesi (kg)" });
    },

    "rocket-performance-form"(data) {
      const chamber = data.result.chamber.mole_fractions.slice(0, 8);
      renderBars(byId("performance-species-chart"), chamber.map(([species, fraction]) => ({
        label: species,
        segments: [{ value: fraction, slot: 1, name: "mol kesri" }],
      })), { unit: "", xLabel: "Mol kesri" });
    },

    "propellant-form"(data) {
      const slots = { liquid: 1, solid: 2, hybrid: 3 };
      const groups = { liquid: "Sıvı–sıvı", solid: "Katı", hybrid: "Hibrit" };
      const rows = [...data.records].sort((a, b) => b.vacuum_isp_max_s - a.vacuum_isp_max_s).map((r) => ({
        label: r.name, min: r.vacuum_isp_min_s, max: r.vacuum_isp_max_s, slot: slots[r.category] || 4, group: groups[r.category] || r.category,
      }));
      renderRanges(byId("propellant-isp-chart"), rows, { xLabel: "Temsilî vakum Isp (s)", unit: "s" });
    },
  };

  document.addEventListener("calculated", (event) => {
    const handler = handlers[event.target.id];
    if (!handler) return;
    Promise.resolve(handler(event.detail, event.target)).catch((error) => {
      console.warn(`Görsel çizilemedi (${event.target.id}):`, error);
    });
  });

  /* ── genel bakış: canlı nozul akış alanı ───────────────────────────── */

  async function overviewHero() {
    const svg = byId("overview-field-chart");
    if (!svg) return;
    try {
      const [geometry, nozzle] = await Promise.all([
        post("/api/v1/nozzle-contour", {
          contour: "rao", throat_area_m2: 0.01, area_ratio: 25, gamma: 1.22, half_angle_deg: 15,
          chamber_radius_ratio: 2.5, bell_length_fraction: 0.8, characteristic_count: 20, throat_angle_deg: 30,
          exit_angle_deg: 8, station_count: 121, axisymmetric_correction: true, include_exports: false,
          include_comparison: false, angular_segments: 32,
        }),
        post("/api/v1/nozzle", {
          chamber_pressure_pa: 7e6, chamber_temperature_k: 3500, throat_area_m2: 0.01, area_ratio: 25,
          ambient_pressure_pa: 0, gamma: 1.22, gas_constant_j_kg_k: 355,
        }),
      ]);
      renderMachField(svg, geometry.stations);
      const set = (id, value) => { const node = byId(id); if (node) node.textContent = value; };
      set("hero-thrust", fmt(nozzle.thrust_n));
      set("hero-isp", fmt(nozzle.specific_impulse_s));
      set("hero-exit-mach", fmt(nozzle.exit_mach));
      set("hero-mdot", fmt(nozzle.mass_flow_kg_s));
    } catch (error) {
      console.warn("Genel bakış görseli çizilemedi:", error);
    }
  }
  overviewHero();

  // transient.js draws its bar, range and histogram charts with the same renderers.
  window.labCharts = { renderBars, renderRanges, renderHistogram };
})();
