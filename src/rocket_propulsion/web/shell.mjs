/* Uygulama kabuğu: görünüm yönlendirme, alt hesap ve grafik sekmeleri, tema,
   hesap motoru durumu ve alan değişince yeniden hesap. Hesap ve çizim app.js'te
   kalır; bu modül yalnız sunuşu düzenler ve formları yeniden gönderir.

   Sunumatik presetleri (vendor/sunumatik, kaynak: PROVENANCE.json):
     · pano-odak.mjs — role=tab sekmelerinin tıklama/klavye davranışı (odakKur)
     · tr-bicim.mjs  — görüntülenen sayılarda Türkçe ondalık virgül
     · palettes.json + palette-library.css — tema seçici                       */

import { odakKur } from "/vendor/sunumatik/pano-odak.mjs";
import "/vendor/sunumatik/tr-bicim.mjs";

const STORAGE = { palette: "rplab.palette", view: "rplab.view" };
const DEFAULT_VIEW = "genel-bakis";
const AUTO_PALETTES = { light: "cerulean-vermilion", dark: "midnight-tangerine" };
// Altı veri serili, datawiz palet_denetimi'nden 0 HATA ile geçen canlı paletler.
const RECOMMENDED = [
  "cerulean-vermilion", "midnight-tangerine",
  "royal-saffron", "heliograph",
  "wada-salvia-indigo", "wada-salvia-indigo-dark",
];

const storage = {
  get(key) { try { return localStorage.getItem(key); } catch { return null; } },
  set(key, value) { try { localStorage.setItem(key, value); } catch { /* özel pencere */ } },
};

let tabSerial = 0;
const nextId = (prefix) => `${prefix}-${++tabSerial}`;

/* ── sekmeler ──────────────────────────────────────────────────────── */

function buildTabList(container, panels, labelOf, extraClass) {
  const list = document.createElement("div");
  list.className = `po-sekmeler ${extraClass}`;
  list.setAttribute("role", "tablist");
  panels.forEach((panel, index) => {
    panel.id ||= nextId("panel");
    panel.setAttribute("role", "tabpanel");
    const tab = document.createElement("button");
    tab.type = "button";
    tab.id = nextId("tab");
    tab.setAttribute("role", "tab");
    tab.setAttribute("aria-controls", panel.id);
    tab.setAttribute("aria-selected", String(index === 0));
    tab.textContent = labelOf(panel);
    panel.setAttribute("aria-labelledby", tab.id);
    list.append(tab);
  });
  container.prepend(list);
}

function buildPaneTabs() {
  for (const set of document.querySelectorAll(".panes")) {
    const panes = [...set.children].filter((child) => child.dataset.pane);
    if (panes.length > 1) buildTabList(set, panes, (pane) => pane.dataset.pane, "pane-tabs");
  }
}

function buildChartTabs() {
  for (const grid of document.querySelectorAll(".chart-grid")) {
    const figures = [...grid.querySelectorAll(":scope > .chart-figure")];
    if (figures.length < 2) continue;
    figures.forEach((figure) => figure.classList.add("po-panel"));
    grid.dataset.tabbed = "true";
    buildTabList(grid, figures,
      (figure) => figure.querySelector("figcaption")?.textContent.trim() || "Grafik", "chart-tabs");
  }
}

/* ── görünüm yönlendirme ───────────────────────────────────────────── */

const views = [...document.querySelectorAll("main > .view")];
const navLinks = [...document.querySelectorAll(".side-nav a[href^='#']")];
const sideNav = document.getElementById("side-nav");
const navToggle = document.querySelector(".nav-toggle");

function viewFor(hash) {
  const id = decodeURIComponent(hash.replace(/^#/, ""));
  if (!id || id === "top") return null;
  const target = document.getElementById(id);
  return target ? target.closest("main > .view") : null;
}

function showView(view, { focus = false } = {}) {
  for (const candidate of views) candidate.hidden = candidate !== view;
  for (const link of navLinks) {
    if (link.getAttribute("href") === `#${view.id}`) link.setAttribute("aria-current", "page");
    else link.removeAttribute("aria-current");
  }
  const heading = view.querySelector("h2");
  document.title = view.id === DEFAULT_VIEW
    ? "Rocket Propulsion Lab"
    : `${heading?.textContent.trim()} · Rocket Propulsion Lab`;
  storage.set(STORAGE.view, view.id);
  window.scrollTo({ top: 0 });
  if (focus && heading) {
    heading.tabIndex = -1;
    heading.focus({ preventScroll: true });
  }
}

function route({ focus = true } = {}) {
  const view = viewFor(location.hash)
    || (!location.hash && viewFor(`#${storage.get(STORAGE.view) || ""}`))
    || document.getElementById(DEFAULT_VIEW);
  showView(view, { focus });
  closeNav();
}

function closeNav() {
  sideNav.dataset.open = "false";
  navToggle.setAttribute("aria-expanded", "false");
}

navToggle.addEventListener("click", () => {
  const open = sideNav.dataset.open !== "true";
  sideNav.dataset.open = String(open);
  navToggle.setAttribute("aria-expanded", String(open));
});

/* ── tema ──────────────────────────────────────────────────────────── */

const themeSelect = document.getElementById("theme-select");
const systemDark = window.matchMedia("(prefers-color-scheme: dark)");
let paletteModes = new Map();

function applyPalette(choice) {
  const id = choice === "auto" || !choice
    ? AUTO_PALETTES[systemDark.matches ? "dark" : "light"]
    : choice;
  const dark = paletteModes.has(id) ? paletteModes.get(id) : id.endsWith("-dark");
  document.documentElement.dataset.palette = id;
  document.documentElement.dataset.mode = dark ? "dark" : "light";
}

function option(value, label) {
  const element = document.createElement("option");
  element.value = value;
  element.textContent = label;
  return element;
}

function group(label, items) {
  const element = document.createElement("optgroup");
  element.label = label;
  element.append(...items);
  return element;
}

async function buildThemePicker() {
  // ?palette=<id> (Sunumatik sayfa kabuğuyla aynı sözleşme) kayıtlı tercihi değiştirmeden uygulanır.
  const stored = new URLSearchParams(location.search).get("palette") || storage.get(STORAGE.palette) || "auto";
  themeSelect.append(option("auto", "Sistemle aynı"));
  try {
    const response = await fetch("/vendor/sunumatik/palettes.json");
    const { palettes } = await response.json();
    paletteModes = new Map(palettes.map((palette) => [palette.id, palette.dark]));
    const byId = new Map(palettes.map((palette) => [palette.id, palette]));
    const label = (palette) => `${palette.ad}${palette.dark && !/koyu/i.test(palette.ad) ? " — koyu" : ""}`;
    const recommended = RECOMMENDED.filter((id) => byId.has(id)).map((id) => option(id, label(byId.get(id))));
    const rest = palettes.filter((palette) => !RECOMMENDED.includes(palette.id))
      .sort((a, b) => a.ad.localeCompare(b.ad, "tr"));
    themeSelect.append(
      group("Önerilen", recommended),
      group("Açık paletler", rest.filter((palette) => !palette.dark).map((palette) => option(palette.id, label(palette)))),
      group("Koyu paletler", rest.filter((palette) => palette.dark).map((palette) => option(palette.id, label(palette)))),
    );
  } catch {
    /* palet listesi yüklenemezse otomatik tema yine çalışır */
  }
  // Liste yüklenemese de istenen palet uygulanır (CSS bütün paletleri zaten taşır).
  if (stored !== "auto" && ![...themeSelect.options].some((item) => item.value === stored)) {
    themeSelect.append(option(stored, stored));
  }
  themeSelect.value = stored;
  applyPalette(themeSelect.value);
}

themeSelect.addEventListener("change", () => {
  storage.set(STORAGE.palette, themeSelect.value);
  applyPalette(themeSelect.value);
});
systemDark.addEventListener("change", () => {
  if (themeSelect.value === "auto") applyPalette("auto");
});

/* ── hesap motoru durumu ───────────────────────────────────────────── */

const status = document.getElementById("engine-status");

async function checkEngine(retry = true) {
  const label = status.querySelector("span");
  try {
    const response = await fetch("/api/v1/health", { cache: "no-store" });
    const body = await response.json();
    if (!response.ok || body.status !== "ok") throw new Error("health");
    status.dataset.state = "ok";
    label.textContent = `Hesap motoru hazır · v${body.version}`;
  } catch {
    // Açılıştaki istek patlamasında tek bir bağlantı reddedilebilir: bir kez daha dene.
    if (retry) { setTimeout(() => checkEngine(false), 1500); return; }
    status.dataset.state = "bad";
    label.textContent = "Sunucuya ulaşılamıyor — app.cmd ile yeniden başlatın";
  }
}

/* ── alan değişince yeniden hesap ──────────────────────────────────── */

const AUTO_FORMS = "form.control-panel, form.catalog-controls";
const timers = new WeakMap();

function resultPanelsOf(form) {
  return [...form.parentElement.querySelectorAll(":scope > .results-panel")];
}

function scheduleRecalculation(form) {
  clearTimeout(timers.get(form));
  timers.set(form, setTimeout(() => {
    const button = form.querySelector("button[type='submit']");
    // Önceki istek sürerken gönderme: yanıtlar sırasız gelip eski sonuç yeniyi ezmesin.
    if (button?.disabled) { scheduleRecalculation(form); return; }
    if (form.checkValidity()) form.requestSubmit();
  }, 250));
}

document.addEventListener("input", (event) => {
  const form = event.target.closest?.(AUTO_FORMS);
  if (form) for (const panel of resultPanelsOf(form)) panel.dataset.stale = "true";
});
document.addEventListener("change", (event) => {
  const form = event.target.closest?.(AUTO_FORMS);
  if (form && event.target.type !== "file") scheduleRecalculation(form);
});

// app.js hesap başlarken form mesajını yalın "form-message" yapar, bitince success/error
// ekler. Ağır hesaplar (tasarım dışı nozul, O/F taraması, kontur) 0,4–2 s sürer: bu
// arada sonuçlar soluk ve aria-busy, bitince yeniden geçerli (Doherty göstergesi).
new MutationObserver((records) => {
  for (const record of records) {
    const message = record.target;
    if (!message.classList.contains("form-message")) continue;
    const form = message.closest(AUTO_FORMS);
    if (!form) continue;
    const done = /\b(success|error)\b/.test(message.className);
    for (const panel of resultPanelsOf(form)) {
      panel.dataset.stale = String(!done);
      panel.setAttribute("aria-busy", String(!done));
    }
  }
}).observe(document.querySelector("main"), { subtree: true, attributes: true, attributeFilter: ["class"] });

/* ── başlat ────────────────────────────────────────────────────────── */

buildPaneTabs();
buildChartTabs();
odakKur();
buildThemePicker();
checkEngine();
setInterval(checkEngine, 30000);
window.addEventListener("hashchange", () => route());
route({ focus: false });
