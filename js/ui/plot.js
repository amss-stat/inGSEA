// ═══════════════════════════════════════════════════════════
//  ui/plot.js  ·  v2.5
//
//  Changes from v2.4:
//  1. Issue (2):  plot-db-link tooltip — already set via title attr in HTML;
//                 updatePlotHeader also sets it programmatically.
//  2. Issue (3):  downloadCurve() — exports current SVG as PNG.
//  3. Issue (8):  drawNESChart() — horizontal bar chart of NES values.
//  4. Issue (9):  exportAllCurves() — batch PNG export of all curves.
//  5. Issue (10): removed Γ/GΓ from renderESStats labels.
//  6. Drag/zoom fixes retained from v2.4.
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

const M       = { l: 64, r: 16, t: 22, b: 48, hitH: 24 };
const SVG_H   = 310;
const MAX_PTS = 4000;

let _zoom = { i0: 0, i1: 1 };
let _cur  = { result: null, pathways: null, geneNames: null, container: null };

// ── Public API ────────────────────────────────────────────────

export function drawCurve(result, pathways, geneNames, container) {
  _cur = { result, pathways, geneNames, container };
  _render();
}

export function resetZoom() {
  _zoom = { i0: 0, i1: 1 };
  _render();
}

export function updatePlotHeader(result) {
  document.getElementById('plot-name').textContent = result.name;
  document.getElementById('plot-section').style.display = 'block';
  const link = document.getElementById('plot-db-link');
  link.href  = msigdbUrl(result.name, result.url);
  // Issue (2): ensure tooltip is set
  link.title = 'Jump to pathway database information';
  link.style.display = 'inline-flex';
  _zoom = { i0: 0, i1: 1 };
}

// ── Issue (3): Download current enrichment curve as PNG ───────
export function downloadCurve(filename) {
  const svgEl = document.getElementById('es-svg');
  if (!svgEl) return;
  _svgToPng(svgEl, filename ?? 'enrichment_curve.png');
}

// ── Issue (8): NES bar chart ──────────────────────────────────
/**
 * Draw a horizontal bar chart of NES values.
 * Positive NES → bar extends right (green).
 * Negative NES → bar extends left (red).
 * Sorted top-to-bottom: largest positive → most negative.
 *
 * @param {Array}  results   — GSEA result objects with .name and .nes
 * @param {Element} container — DOM element to render into
 */
export function drawNESChart(results, container) {
  if (!results?.length || !container) return;

  // Sort: largest NES at top
  const sorted = [...results].sort((a, b) => b.nes - a.nes);
  const n      = sorted.length;

  const BAR_H  = Math.max(14, Math.min(28, Math.floor(520 / n)));
  const GAP    = 3;
  const ML     = 180;   // left margin for pathway names
  const MR     = 50;    // right margin for value labels
  const MT     = 28;    // top margin
  const MB     = 24;    // bottom margin
  const W      = Math.max(container.clientWidth || 680, 500);
  const plotW  = W - ML - MR;
  const H      = MT + MB + n * (BAR_H + GAP);

  // Axis scale: symmetric around 0
  const absMax = Math.max(...sorted.map(r => Math.abs(r.nes)), 0.1);
  const xScale = (plotW / 2) / absMax;   // px per NES unit
  const x0     = ML + plotW / 2;         // x-coordinate of NES = 0

  const toX = nes => x0 + nes * xScale;

  // Axis ticks (symmetric, 5 steps each side)
  const tickStep = _niceStep(absMax / 4);
  const ticks    = [];
  for (let v = 0; v <= absMax + tickStep * 0.5; v += tickStep)
    ticks.push(v, -v);
  const uniqueTicks = [...new Set(ticks)].sort((a, b) => a - b);

  const tickSVG = uniqueTicks.map(v => {
    const x = _f(toX(v));
    return `
    <line x1="$${x}" y1="$${MT - 6}" x2="$${x}" y2="$${MT + n*(BAR_H+GAP)}"
          stroke="#e4e8ee" stroke-width="0.8" stroke-dasharray="3,4"/>
    <text x="$${x}" y="$${MT - 8}" text-anchor="middle"
          fill="#8090a8" font-size="9">${v.toFixed(1)}</text>`;
  }).join('');

  // Bars
  const barsSVG = sorted.map((r, i) => {
    const y    = MT + i * (BAR_H + GAP);
    const nes  = r.nes;
    const barX = nes >= 0 ? x0 : toX(nes);
    const barW = Math.abs(nes) * xScale;
    const col  = nes >= 0 ? '#1c6e41' : '#b01c1c';
    const lCol = nes >= 0 ? '#134d2e' : '#7a1313';
    const sig  = r.pCauchy != null && r.pCauchy < 0.05;
    const nameX = ML - 5;

    // Significance marker
    const sigMark = sig ? '★' : '';

    // Clip name to fit in left margin
    const shortName = _esc(r.name.length > 28
      ? r.name.slice(0, 26) + '…'
      : r.name);

    return `
    <g class="nes-bar-group" data-name="${_esc(r.name)}">
      <rect x="$${_f(barX)}" y="$${y}" width="$${_f(Math.max(barW, 1))}" height="$${BAR_H}"
            fill="${col}" opacity="0.82" rx="2"/>
      <text x="$${nameX}" y="$${y + BAR_H/2 + 4}"
            text-anchor="end" fill="${sig ? '#1a3a6a' : '#4a5568'}"
            font-size="${Math.min(11, BAR_H - 3)}"
            font-weight="$${sig ? '600' : '400'}">$${shortName} ${sigMark}</text>
      <text x="$${_f(toX(nes) + (nes >= 0 ? 3 : -3))}" y="$${y + BAR_H/2 + 4}"
            text-anchor="${nes >= 0 ? 'start' : 'end'}"
            fill="${lCol}" font-size="9" font-family="'JetBrains Mono',monospace"
            font-weight="500">${nes.toFixed(2)}</text>
    </g>`;
  }).join('');

  // Zero line
  const zeroSVG = `
  <line x1="$${_f(x0)}" y1="$${MT - 10}"
        x2="$${_f(x0)}" y2="$${MT + n*(BAR_H+GAP)}"
        stroke="#6b7280" stroke-width="1.5"/>`;

  // Axis label
  const axisLabel = `
  <text x="$${_f(x0)}" y="$${H - 4}" text-anchor="middle"
        fill="#7a8698" font-size="11" font-weight="500">
    Normalized Enrichment Score (NES)
  </text>
  <text x="$${_f(x0 - plotW/4)}" y="$${H - 4}" text-anchor="middle"
        fill="#b01c1c" font-size="9.5">◀ Negative</text>
  <text x="$${_f(x0 + plotW/4)}" y="$${H - 4}" text-anchor="middle"
        fill="#1c6e41" font-size="9.5">Positive ▶</text>
  <text x="$${W - MR + 2}" y="$${MT + 10}" text-anchor="end"
        fill="#9aa3b0" font-size="8.5" font-style="italic">
    ★ p<sub>Cauchy</sub>&lt;0.05
  </text>`;

  const svg = `<svg xmlns="http://www.w3.org/2000/svg"
  viewBox="0 0 $${W} $${H}" width="$${W}" height="$${H}"
  id="nes-svg" style="font-family:'Inter',sans-serif"
  role="img" aria-label="NES bar chart">
  <rect width="$${W}" height="$${H}" fill="#fafbfc"/>
  ${tickSVG}
  ${zeroSVG}
  ${barsSVG}
  ${axisLabel}
</svg>`;

  container.innerHTML = svg;
}

// ── Issue (9): Batch export all enrichment curves ─────────────
/**
 * Export all pathway enrichment curves as individual PNG files,
 * packaged into a ZIP using the JSZip library if available,
 * or downloaded sequentially (one per second) if not.
 *
 * @param {Array}   results   — all GSEA results
 * @param {Array}   pathways  — pathway list (with masks)
 * @param {Array}   geneNames
 * @param {Function} onProgress  — callback(done, total)
 */
export async function exportAllCurves(results, pathways, geneNames, onProgress) {
  if (!results?.length) return;

  // Use an off-screen container sized to a standard width
  const offscreen = document.createElement('div');
  offscreen.style.cssText =
    'position:fixed;left:-9999px;top:0;width:720px;height:400px;' +
    'background:white;visibility:hidden;';
  document.body.appendChild(offscreen);

  const total = results.length;
  const pngs  = [];   // [{ name, blob }]

  for (let i = 0; i < total; i++) {
    const r = results[i];

    // Save and restore zoom
    const savedZoom = { ..._zoom
