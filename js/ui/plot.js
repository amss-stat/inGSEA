// ═══════════════════════════════════════════════════════════
//  ui/plot.js  ·  v2.5
//  Changes:
//  • renderESStats: removed Γ/GΓ → plain "par" label
//  • Added exportCurvePNG(svgEl, filename) for single download
//  • Added exportAllCurves(results, pathways, geneNames) for batch
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

const M       = { l: 64, r: 16, t: 22, b: 48, hitH: 24 };
const SVG_H   = 310;
const MAX_PTS = 4000;

let _zoom = { i0: 0, i1: 1 };
let _cur  = { result: null, pathways: null, geneNames: null, container: null };

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
  link.style.display = 'inline-flex';
  _zoom = { i0: 0, i1: 1 };
}

// ── Single-curve PNG download (issue 3) ───────────────────────
/**
 * Export the currently displayed SVG as a PNG file.
 * @param {string} filename
 */
export function exportCurrentCurvePNG(filename) {
  const svgEl = document.querySelector('#svg-wrap #es-svg');
  if (!svgEl) return;
  _svgToPNG(svgEl, filename || 'enrichment_curve.png');
}

/**
 * Batch-export all results as PNG files (issue 9).
 * Renders each curve off-screen, downloads sequentially with
 * a small delay to avoid browser download blocking.
 *
 * @param {Array}  results
 * @param {Array}  pathways
 * @param {Array}  geneNames
 * @param {number} [width=900]
 */
export async function exportAllCurves(results, pathways, geneNames, width = 900) {
  const offscreen = document.createElement('div');
  offscreen.style.cssText =
    'position:fixed;left:-9999px;top:0;width:' + width + 'px;visibility:hidden';
  document.body.appendChild(offscreen);

  for (let i = 0; i < results.length; i++) {
    const r = results[i];
    // Render the SVG into the offscreen div
    _renderInto(r, pathways, geneNames, offscreen, width);
    const svgEl = offscreen.querySelector('#es-svg');
    if (!svgEl) continue;
    const safe = r.name.replace(/[^a-z0-9_\-]/gi, '_').slice(0, 80);
    await _svgToPNGAsync(svgEl, `igsea_${String(i + 1).padStart(3, '0')}_${safe}.png`);
    // Small pause between downloads
    await _sleep(120);
  }

  document.body.removeChild(offscreen);
}

// ── Internal renderer ─────────────────────────────────────────
function _render() {
  const { result, pathways, geneNames, container } = _cur;
  if (!result || !container) return;
  const W = Math.max(container.clientWidth || 680, 400);
  _renderInto(result, pathways, geneNames, container, W);
  _attachInteraction(
    container,
    result.curve,
    result.obsOrd,
    geneNames,
    _buildLayout(W, result.curve.length)
  );
}

/** Build layout constants for a given width and gene count. */
function _buildLayout(W, nG) {
  const pw   = W - M.l - M.r;
  const ph   = SVG_H - M.t - M.b - M.hitH - 6;
  const r0   = Math.max(0,    Math.floor(_zoom.i0 * (nG - 1)));
  const r1   = Math.min(nG-1, Math.ceil (_zoom.i1 * (nG - 1)));
  const nVis = r1 - r0 + 1;
  const toX  = r  => M.l + ((r - r0) / Math.max(nVis - 1, 1)) * pw;
  return { W, pw, ph, r0, r1, nVis, toX };
}

function _renderInto(result, pathways, geneNames, container, W) {
  W = Math.max(W || 680, 400);
  const pw   = W - M.l - M.r;
  const ph   = SVG_H - M.t - M.b - M.hitH - 6;

  const curve  = result.curve;
  const nG     = curve.length;
  const ord    = result.obsOrd;
  const isPos  = result.es >= 0;
  const lineC  = isPos ? '#1c6e41' : '#b01c1c';
  const fillC  = isPos ? 'rgba(28,110,65,.11)' : 'rgba(176,28,28,.09)';

  const r0   = Math.max(0,    Math.floor(_zoom.i0 * (nG - 1)));
  const r1   = Math.min(nG-1, Math.ceil (_zoom.i1 * (nG - 1)));
  const nVis = r1 - r0 + 1;

  let lo = 0, hi = 0;
  for (let i = r0; i <= r1; i++) {
    if (curve[i] < lo) lo = curve[i];
    if (curve[i] > hi) hi = curve[i];
  }
  const pad  = Math.max(-lo, hi) * 0.15 + 0.02;
  const yMin = lo - pad, yMax = hi + pad, yR = yMax - yMin;

  const toX = r  => M.l + ((r - r0) / Math.max(nVis - 1, 1)) * pw;
  const toY = es => M.t + ph * (1 - (es - yMin) / yR);
  const y0  = toY(0);

  const pkIdx = result.peakIdx;
  const step  = Math.max(1, Math.ceil(nVis / MAX_PTS));

  let pathD = `M${_f(toX(r0))},${_f(toY(curve[r0]))}`;
  let fillD = `M${_f(toX(r0))},${_f(y0)}`;

  for (let i = r0 + step; i <= r1; i += step) {
    if (pkIdx >= r0 && pkIdx <= r1 && i > pkIdx && i - step <= pkIdx) {
      pathD += ` L${_f(toX(pkIdx))},${_f(toY(curve[pkIdx]))}`;
      fillD += ` L${_f(toX(pkIdx))},${_f(toY(curve[pkIdx]))}`;
    }
    pathD += ` L${_f(toX(i))},${_f(toY(curve[i]))}`;
    fillD += ` L${_f(toX(i))},${_f(toY(curve[i]))}`;
  }
  pathD += ` L${_f(toX(r1))},${_f(toY(curve[r1]))}`;
  fillD += ` L${_f(toX(r1))},${_f(toY(curve[r1]))} L${_f(toX(r1))},${_f(y0)} Z`;

  const ticks   = Array.from({ length: 6 }, (_, k) => yMin + yR * k / 5);
  const tickSVG = ticks.map(v => {
    const y = _f(toY(v));
    return `<line x1="${M.l}" y1="${y}" x2="${M.l+pw}" y2="${y}"
      stroke="#e4e8ee" stroke-width="0.7" stroke-dasharray="3,4"/>
    <text x="${M.l-5}" y="${_f(toY(v)+3.5)}"
      text-anchor="end" fill="#8090a8" font-size="9.5">${v.toFixed(2)}</text>`;
  }).join('');

  const pathway = pathways.find(p => p.name === result.name);
  const mask    = pathway?.mask ?? null;
  const hitY    = M.t + ph + 6;
  let   hitSVG  = '';
  if (mask) {
    const parts = [];
    for (let i = r0; i <= r1; i++)
      if (mask[ord[i]]) parts.push(`M${_f(toX(i))},${hitY} v${M.hitH-5}`);
    if (parts.length)
      hitSVG = `<path d="${parts.join(' ')}"
        stroke="${lineC}" stroke-width="1.1" opacity="0.65" fill="none"/>`;
  }

  let peakSVG = '';
  if (pkIdx >= r0 && pkIdx <= r1) {
    const px = _f(toX(pkIdx)), py = _f(toY(curve[pkIdx]));
    peakSVG = `
    <line x1="${px}" y1="${M.t}" x2="${px}" y2="${M.t+ph}"
      stroke="#a05c07" stroke-width="1" stroke-dasharray="4,3"
      opacity="0.7" clip-path="url(#cp)"/>
    <circle cx="${px}" cy="${py}" r="4"
      fill="#a05c07" stroke="white" stroke-width="1.2"
      clip-path="url(#cp)"/>`;
  }

  const xMid     = _f(M.l + pw / 2);
  const isZoomed = _zoom.i0 > 0.001 || _zoom.i1 < 0.999;
  const zoomTxt = isZoomed
    ? `Rank ${r0 + 1} – ${r1 + 1}`
    : '';

  const svg = `<svg xmlns="http://www.w3.org/2000/svg"
  viewBox="0 0 ${W} ${SVG_H}" width="${W}" height="${SVG_H}"
  id="es-svg" style="font-family:'Inter',sans-serif;overflow:visible"
  role="img" aria-label="Enrichment plot: ${_esc(result.name)}">

  <defs>
    <clipPath id="cp">
      <rect x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"/>
    </clipPath>
  </defs>

  <rect width="${W}" height="${SVG_H}" fill="#fafbfc"/>
  <rect x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="white" stroke="#dde2ea" stroke-width="0.8"/>

  ${tickSVG}
  ${[1,2,3,4].map(k => {
    const x = _f(M.l + pw * k / 5);
    return `<line x1="${x}" y1="${M.t}" x2="${x}" y2="${M.t+ph}"
      stroke="#e4e8ee" stroke-width="0.7" stroke-dasharray="3,4"/>`;
  }).join('')}

  <line x1="${M.l}" y1="${_f(y0)}" x2="${M.l+pw}" y2="${_f(y0)}"
        stroke="#bbc4d0" stroke-width="0.9"/>

  <path d="${fillD}" fill="${fillC}" clip-path="url(#cp)"/>
  <path id="es-path" d="${pathD}" fill="none" stroke="${lineC}"
        stroke-width="1.9" stroke-linejoin="round" stroke-linecap="round"
        clip-path="url(#cp)"/>

  ${peakSVG}

  <rect x="${M.l}" y="${hitY}" width="${pw}" height="${M.hitH}"
        fill="#f0f2f6" stroke="#dde2ea" stroke-width="0.6"/>
  ${hitSVG}
  <text x="${xMid}" y="${hitY+M.hitH-5}" text-anchor="middle"
        fill="#9aa3b0" font-size="8" font-weight="600"
        letter-spacing=".06em">GENE HITS</text>

  <line x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t+ph}"
        stroke="#adb8c4" stroke-width="1.3"/>
  <line x1="${M.l}" y1="${M.t+ph}" x2="${M.l+pw}" y2="${M.t+ph}"
        stroke="#adb8c4" stroke-width="1.3"/>

  <text transform="translate(13,${_f(M.t+ph/2)}) rotate(-90)"
        text-anchor="middle" fill="#7a8698"
        font-size="10.5" font-weight="500">Enrichment Score</text>

  <text x="${_f(toX(r0))}" y="${hitY+M.hitH+14}"
        text-anchor="start" fill="#7a8698" font-size="9.5">${r0+1}</text>
  <text x="${xMid}" y="${hitY+M.hitH+14}"
        text-anchor="middle" fill="#7a8698" font-size="9.5">
    ${Math.round((r0+r1)/2)+1}
  </text>
  <text x="${_f(toX(r1))}" y="${hitY+M.hitH+14}"
        text-anchor="end" fill="#7a8698" font-size="9.5">${r1+1}</text>
  <text x="${xMid}" y="${hitY+M.hitH+30}"
        text-anchor="middle" fill="#7a8698"
        font-size="11" font-weight="500">Gene rank</text>

  <text x="${M.l+pw-2}" y="${M.t+13}" text-anchor="end"
        fill="${isZoomed ? '#a05c07' : '#9aa3b0'}"
        font-size="8.5" font-style="italic">${zoomTxt}</text>

  <rect id="hover-overlay"
        x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="transparent" style="cursor:crosshair"/>

  <line id="h-line"
        x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t+ph}"
        stroke="#1a5fa8" stroke-width="1" stroke-dasharray="3,3"
        opacity="0" pointer-events="none"/>
  <circle id="h-dot" cx="${M.l}" cy="${M.t}" r="3.5"
          fill="#1a5fa8" stroke="white" stroke-width="1.2"
          opacity="0" pointer-events="none"/>
</svg>`;

  container.innerHTML = svg;
}

// ── Interaction ───────────────────────────────────────────────
function _attachInteraction(container, curve, ord, geneNames,
                             { pw, ph, r0, r1, nVis, toX }) {
  // Re-derive toY from current zoom state
  let lo = 0, hi = 0;
  for (let i = r0; i <= r1; i++) {
    if (curve[i] < lo) lo = curve[i];
    if (curve[i] > hi) hi = curve[i];
  }
  const pad  = Math.max(-lo, hi) * 0.15 + 0.02;
  const yMin = lo - pad, yMax = hi + pad, yR = yMax - yMin;
  const toY  = es => M.t + (SVG_H - M.t - M.b - M.hitH - 6) * (1 - (es - yMin) / yR);

  const svgEl   = container.querySelector('#es-svg');
  const overlay = container.querySelector('#hover-overlay');
  const hLine   = container.querySelector('#h-line');
  const hDot    = container.querySelector('#h-dot');
  const tooltip = document.getElementById('svg-tooltip');
  if (!svgEl || !overlay || !tooltip) return;

  let _ctm = null;
  const getCTM  = () => { if (!_ctm) _ctm = svgEl.getScreenCTM(); return _ctm; };
  const toSVGPt = (cx, cy) =>
    new DOMPoint(cx, cy).matrixTransform(getCTM().inverse());

  const svgXtoRank = x =>
    Math.max(r0, Math.min(r1,
      r0 + Math.round((x - M.l) / pw * (nVis - 1))
    ));

  let _isDragging = false;

  overlay.addEventListener('mousemove', e => {
    if (_isDragging) return;
    const pt   = toSVGPt(e.clientX, e.clientY);
    const rank = svgXtoRank(pt.x);
    const es   = curve[rank];
    const lx   = _f(toX(rank));
    const ly   = _f(toY(es));

    hLine.setAttribute('x1', lx); hLine.setAttribute('x2', lx);
    hLine.setAttribute('opacity', '0.55');
    hDot.setAttribute('cx', lx);  hDot.setAttribute('cy', ly);
    hDot.setAttribute('opacity', '1');

    tooltip.style.display = 'block';
    tooltip.style.left    = (e.clientX + 14) + 'px';
    tooltip.style.top     = (e.clientY - 42) + 'px';
    tooltip.innerHTML =
      `<strong>${_esc(geneNames[ord[rank]])}</strong><br>` +
      `Rank #${rank+1}<br>ES&nbsp;${es.toFixed(4)}`;
  });

  overlay.addEventListener('mouseleave', () => {
    if (_isDragging) return;
    _hideHover();
  });

  const _hideHover = () => {
    hLine.setAttribute('opacity', '0');
    hDot.setAttribute('opacity', '0');
    tooltip.style.display = 'none';
  };

  overlay.addEventListener('wheel', e => {
    e.preventDefault();
    _ctm = null;
    const pt      = toSVGPt(e.clientX, e.clientY);
    const foc     = Math.max(0, Math.min(1, (pt.x - M.l) / pw));
    const zFactor = e.deltaY < 0 ? 0.7 : 1 / 0.7;
    const span    = _zoom.i1 - _zoom.i0;
    const newSpan = Math.max(0.01, Math.min(1, span * zFactor));
    let   i0 = _zoom.i0 + foc * (span - newSpan);
    let   i1 = i0 + newSpan;
    if (i0 < 0) { i0 = 0; i1 = newSpan; }
    if (i1 > 1) { i1 = 1; i0 = 1 - newSpan; }
    _zoom = { i0, i1 };
    _render();
  }, { passive: false });

  let _drag = null;

  overlay.addEventListener('mousedown', e => {
    if (e.button !== 0) return;
    e.preventDefault();
    _ctm = null;
    const pt = toSVGPt(e.clientX, e.clientY);
    _drag = { startSVGX: pt.x, startI0: _zoom.i0, startI1: _zoom.i1 };
    _isDragging = false;
    overlay.style.cursor = 'grab';
  });

  const _onMove = e => {
    if (!_drag) return;
    _ctm = null;
    const pt   = toSVGPt(e.clientX, e.clientY);
    const dSVG = pt.x - _drag.startSVGX;

    if (Math.abs(dSVG) > 3) {
      _isDragging = true;
      _hideHover();
      overlay.style.cursor = 'grabbing';
    }
    if (!_isDragging) return;

    const span   = _drag.startI1 - _drag.startI0;
    const dFrac  = -(dSVG / pw) * span;
    let   i0 = _drag.startI0 + dFrac;
    let   i1 = _drag.startI1 + dFrac;
    if (i0 < 0) { i0 = 0; i1 = span; }
    if (i1 > 1) { i1 = 1; i0 = 1 - span; }
    _zoom = { i0, i1 };
  };

  const _onUp = () => {
    if (!_drag) return;
    const wasDragging = _isDragging;
    _drag = null; _isDragging = false;
    overlay.style.cursor = 'crosshair';
    if (wasDragging) _render();
  };

  window.addEventListener('mousemove', _onMove, { passive: true });
  window.addEventListener('mouseup',   _onUp);

  const obs = new MutationObserver(() => {
    window.removeEventListener('mousemove', _onMove);
    window.removeEventListener('mouseup',   _onUp);
    obs.disconnect();
  });
  obs.observe(container, { childList: true });

  overlay.addEventListener('dblclick', () => {
    _zoom = { i0: 0, i1: 1 };
    _ctm  = null;
    _render();
  });

  const _onResize = () => { _ctm = null; };
  window.addEventListener('resize', _onResize, { passive: true });
  const obs2 = new MutationObserver(() => {
    window.removeEventListener('resize', _onResize);
    obs2.disconnect();
  });
  obs2.observe(container, { childList: true });
}

// ── ES stats row ──────────────────────────────────────────────
export function renderESStats(result, engine, showFDR) {
  const fmt  = p => p == null ? '—'
    : Math.abs(p) < 0.001 ? p.toExponential(2) : p.toFixed(4);
  const sc   = p => p != null && p < 0.05 ? 'sig' : '';
  const isPar = engine === 'parametric';

  // Issue (10): no Γ/GΓ — use plain "(par)" label
  const cells = [
    { l: 'ES',      v: result.es.toFixed(4),    c: result.es    >= 0 ? 'pos' : 'neg' },
    { l: 'NES',     v: result.nes.toFixed(3),   c: result.nes   >= 0 ? 'pos' : 'neg' },
    { l: 'NES-AD',  v: result.nes_ad.toFixed(3), c: '' },
    { l: isPar ? 'p<sub>KS</sub>&thinsp;(par)'  : 'p<sub>KS</sub>',
      v: fmt(result.pKS),     c: sc(result.pKS) },
    { l: isPar ? 'p<sub>AD</sub>&thinsp;(par)'  : 'p<sub>AD</sub>',
      v: fmt(result.pAD),     c: sc(result.pAD) },
    { l: 'p<sub>Cauchy</sub>',
      v: fmt(result.pCauchy), c: sc(result.pCauchy) },
    ...(showFDR && result.fdr != null
      ? [{ l: 'FDR (BH)', v: fmt(result.fdr), c: sc(result.fdr) }]
      : []),
    { l: 'Size', v: `${result.size}`, c: '' }
  ];

  document.getElementById('es-stats').innerHTML = cells.map(c =>
    `<div class="es-cell">
       <span class="es-lbl">${c.l}</span>
       <span class="es-val ${c.c}">${c.v}</span>
     </div>`
  ).join('');
}

// ── NES bar chart (issue 8) ───────────────────────────────────
/**
 * Draw a horizontal NES bar chart for all results.
 * Positive NES → bars face right (green).
 * Negative NES → bars face left (red).
 * Sorted: positive (largest first) on top, negative below.
 *
 * @param {Array}       results   array of result objects with .name and .nes
 * @param {HTMLElement} container
 */
export function drawNESChart(results, container) {
  if (!results || results.length === 0) return;

  const sorted = [...results].sort((a, b) => b.nes - a.nes);

  // Issue (3): layout constants tuned for long names
  const barH   = 20;
  const gap    = 4;
  const padT   = 20;
  const padB   = 32;
  const padL   = 8;
  const padR   = 12;
  // Value label width (right of bar)
  const valW   = 52;

  const nBars  = sorted.length;
  const totalH = padT + nBars * (barH + gap) + padB;

  const W      = Math.max(container.clientWidth || 700, 500);

  // Issue (3): label area is 40% of total width, capped
  const labelW = Math.min(Math.floor(W * 0.40), 260);
  const chartW = W - padL - labelW - valW - padR;

  const maxAbs = sorted.reduce((m, r) => Math.max(m, Math.abs(r.nes)), 0) || 1;
  const zeroX  = padL + labelW + chartW / 2;
  const scale  = v => (v / maxAbs) * (chartW / 2);

  const barRows = sorted.map((r, i) => {
    const y      = padT + i * (barH + gap);
    const yMid   = y + barH / 2 + 4;
    const nes    = r.nes;
    const isPos  = nes >= 0;
    const barPx  = Math.max(scale(Math.abs(nes)), 1);
    const barX   = isPos ? zeroX : zeroX - barPx;
    const fill   = isPos ? '#1c6e41' : '#b01c1c';
    const fillBg = isPos ? 'rgba(28,110,65,.08)' : 'rgba(176,28,28,.07)';

    const star = r.pCauchy < 0.001 ? '***'
               : r.pCauchy < 0.01  ? '**'
               : r.pCauchy < 0.05  ? '*' : '';

    // Issue (3): label rendered as SVG foreignObject for CSS ellipsis,
    // falling back to SVG text with title for tooltip
    const labelX   = padL + 4;
    const labelMaxW = labelW - 10;

    // Value shown to the right (positive) or left (negative) of bar
    const nesStr   = nes.toFixed(2) + (star ? ' ' + star : '');
    const valX     = isPos
      ? zeroX + barPx + 4
      : zeroX - barPx - 4;
    const valAnchor = isPos ? 'start' : 'end';

    return `
      <rect x="${padL}" y="${y}" width="${labelW - 4}" height="${barH}"
            fill="${fillBg}" rx="2"/>
      <text x="${_f(labelX + labelMaxW)}" y="${_f(yMid)}"
            text-anchor="end" font-size="10.5"
            fill="#2a3040" font-family="'Inter',sans-serif">
        <title>${_esc(r.name)}</title>
        ${_esc(_truncate(r.name, 34))}
      </text>
      <rect x="${_f(barX)}" y="${_f(y + 2)}"
            width="${_f(barPx)}" height="${barH - 4}"
            fill="${fill}" opacity="0.80" rx="2"/>
      <text x="${_f(valX)}" y="${_f(yMid)}"
            text-anchor="${valAnchor}"
            font-size="9.5" font-weight="500"
            fill="${fill}"
            font-family="'JetBrains Mono',monospace">
        ${_esc(nesStr)}
      </text>`;
  }).join('');

  // Issue (3): zero line only, no tick marks
  const zeroLine = `
    <line x1="${_f(zeroX)}" y1="${padT - 4}"
          x2="${_f(zeroX)}" y2="${padT + nBars*(barH+gap)}"
          stroke="#adb8c4" stroke-width="1.4"/>`;

  // Legend
  const legendY = padT + nBars * (barH + gap) + 18;
  const legend  = `
    <text x="${_f(W / 2)}" y="${legendY}"
          text-anchor="middle" font-size="9" fill="#8090a8"
          font-family="'Inter',sans-serif">
      NES · * p&lt;0.05 · ** p&lt;0.01 · *** p&lt;0.001 (Cauchy)
    </text>`;

  const svg = `<svg xmlns="http://www.w3.org/2000/svg"
  viewBox="0 0 ${W} ${totalH}" width="${W}" height="${totalH}"
  id="nes-svg" style="font-family:'Inter',sans-serif"
  role="img" aria-label="NES bar chart">
  <rect width="${W}" height="${totalH}" fill="#fafbfc"/>
  ${zeroLine}
  ${barRows}
  ${legend}
</svg>`;

  container.innerHTML = svg;
}

// Truncate string to maxLen characters with ellipsis
function _truncate(s, maxLen) {
  if (!s) return '';
  return s.length <= maxLen ? s : s.slice(0, maxLen - 1) + '…';
}

// ── PNG export helpers ────────────────────────────────────────
/**
 * Convert an SVG element to a PNG and trigger download.
 * Embeds Google Fonts via a style tag (best-effort; fallback to system fonts).
 */
function _svgToPNG(svgEl, filename) {
  _svgToPNGAsync(svgEl, filename).catch(console.error);
}

function _svgToPNGAsync(svgEl, filename) {
  return new Promise((resolve, reject) => {
    const svgClone = svgEl.cloneNode(true);
    // Embed a minimal font declaration so the PNG is self-contained
    const style = document.createElementNS('http://www.w3.org/2000/svg', 'style');
    style.textContent = `
      text { font-family: 'Inter', Arial, sans-serif; }
      text[font-family*="Mono"], text[font-family*="mono"] {
        font-family: 'Courier New', monospace;
      }`;
    svgClone.insertBefore(style, svgClone.firstChild);

    const W   = parseInt(svgEl.getAttribute('width'),  10) || 800;
    const H   = parseInt(svgEl.getAttribute('height'), 10) || 400;
    const scale = 2;   // 2× for retina / print quality
    svgClone.setAttribute('width',  W * scale);
    svgClone.setAttribute('height', H * scale);

    const blob = new Blob(
      [new XMLSerializer().serializeToString(svgClone)],
      { type: 'image/svg+xml' }
    );
    const url = URL.createObjectURL(blob);
    const img = new Image();
    img.onload = () => {
      const canvas  = document.createElement('canvas');
      canvas.width  = W * scale;
      canvas.height = H * scale;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = '#ffffff';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.drawImage(img, 0, 0);
      URL.revokeObjectURL(url);
      canvas.toBlob(pngBlob => {
        if (!pngBlob) { reject(new Error('Canvas toBlob failed')); return; }
        const a = Object.assign(document.createElement('a'), {
          href:     URL.createObjectURL(pngBlob),
          download: filename
        });
        a.click();
        setTimeout(() => URL.revokeObjectURL(a.href), 5000);
        resolve();
      }, 'image/png');
    };
    img.onerror = () => { URL.revokeObjectURL(url); reject(new Error('SVG load failed')); };
    img.src = url;
  });
}

const _sleep = ms => new Promise(r => setTimeout(r, ms));

const _f   = v => v.toFixed(2);
const _esc = s => (s ?? '')
  .replace(/&/g, '&amp;')
  .replace(/</g, '&lt;')
  .replace(/>/g, '&gt;');
