// ═══════════════════════════════════════════════════════════
//  ui/plot.js  ·  v2.4-final
//
//  Fixes:
//  1. Drag re-render removed listener bug:
//     Pan state is applied immediately to _zoom but re-render
//     is deferred to mouseup (or throttled with rAF).
//  2. Removed duplicate id attribute on hover-overlay.
//  3. CTM is re-fetched after each re-render (new SVG element).
//  4. Drag direction: moving right reveals earlier ranks (correct).
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

// ── Renderer ──────────────────────────────────────────────────
function _render() {
  const { result, pathways, geneNames, container } = _cur;
  if (!result || !container) return;

  const W  = Math.max(container.clientWidth || 680, 400);
  const pw = W - M.l - M.r;
  const ph = SVG_H - M.t - M.b - M.hitH - 6;

  const curve = result.curve;
  const nG    = curve.length;
  const ord   = result.obsOrd;
  const isPos = result.es >= 0;
  const lineC = isPos ? '#1c6e41' : '#b01c1c';
  const fillC = isPos ? 'rgba(28,110,65,.11)' : 'rgba(176,28,28,.09)';

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

  // All coordinate functions work in SVG viewBox units
  const toX = r  => M.l + ((r - r0) / Math.max(nVis - 1, 1)) * pw;
  const toY = es => M.t + ph * (1 - (es - yMin) / yR);
  const y0  = toY(0);

  // Build path with subsampling (always include peak)
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
  // Always draw last point
  pathD += ` L${_f(toX(r1))},${_f(toY(curve[r1]))}`;
  fillD += ` L${_f(toX(r1))},${_f(toY(curve[r1]))} L${_f(toX(r1))},${_f(y0)} Z`;

  // Y axis ticks
  const ticks  = Array.from({ length: 6 }, (_, k) => yMin + yR * k / 5);
  const tickSVG = ticks.map(v => {
    const y = _f(toY(v));
    return `<line x1="${M.l}" y1="${y}" x2="${M.l+pw}" y2="${y}"
      stroke="#e4e8ee" stroke-width="0.7" stroke-dasharray="3,4"/>
    <text x="${M.l-5}" y="${_f(toY(v)+3.5)}"
      text-anchor="end" fill="#8090a8" font-size="9.5">${v.toFixed(2)}</text>`;
  }).join('');

  // Gene-hit strips (visible window only)
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

  // Peak marker (only if in visible window)
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
  const zoomTxt  = isZoomed
    ? `Rank ${r0+1}–${r1+1} · drag to pan · scroll to zoom · dbl-click reset`
    : `Scroll to zoom · drag to pan · dbl-click reset`;

  const svg = `<svg xmlns="http://www.w3.org/2000/svg"
  viewBox="0 0 ${W} ${SVG_H}" width="${W}" height="${SVG_H}"
  id="es-svg" style="font-family:'Inter',sans-serif;overflow:visible"
  role="img" aria-label="Enrichment plot: ${_esc(result.name)}">

  <defs>
    <clipPath id="cp">
      <rect x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"/>
    </clipPath>
  </defs>

  <!-- Background -->
  <rect width="${W}" height="${SVG_H}" fill="#fafbfc"/>
  <rect x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="white" stroke="#dde2ea" stroke-width="0.8"/>

  <!-- Grid -->
  ${tickSVG}
  ${[1,2,3,4].map(k=>{const x=_f(M.l+pw*k/5);
    return `<line x1="${x}" y1="${M.t}" x2="${x}" y2="${M.t+ph}"
      stroke="#e4e8ee" stroke-width="0.7" stroke-dasharray="3,4"/>`;
  }).join('')}

  <!-- Zero line -->
  <line x1="${M.l}" y1="${_f(y0)}" x2="${M.l+pw}" y2="${_f(y0)}"
        stroke="#bbc4d0" stroke-width="0.9"/>

  <!-- Fill & curve -->
  <path d="${fillD}" fill="${fillC}" clip-path="url(#cp)"/>
  <path id="es-path" d="${pathD}" fill="none" stroke="${lineC}"
        stroke-width="1.9" stroke-linejoin="round" stroke-linecap="round"
        clip-path="url(#cp)"/>

  <!-- Peak -->
  ${peakSVG}

  <!-- Gene-hit strip -->
  <rect x="${M.l}" y="${hitY}" width="${pw}" height="${M.hitH}"
        fill="#f0f2f6" stroke="#dde2ea" stroke-width="0.6"/>
  ${hitSVG}
  <text x="${xMid}" y="${hitY+M.hitH-5}" text-anchor="middle"
        fill="#9aa3b0" font-size="8" font-weight="600"
        letter-spacing=".06em">GENE HITS</text>

  <!-- Axes -->
  <line x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t+ph}"
        stroke="#adb8c4" stroke-width="1.3"/>
  <line x1="${M.l}" y1="${M.t+ph}" x2="${M.l+pw}" y2="${M.t+ph}"
        stroke="#adb8c4" stroke-width="1.3"/>

  <!-- Axis titles -->
  <text transform="translate(13,${_f(M.t+ph/2)}) rotate(-90)"
        text-anchor="middle" fill="#7a8698"
        font-size="10.5" font-weight="500">Enrichment Score</text>

  <!-- X-axis labels -->
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

  <!-- Zoom status hint -->
  <text x="${M.l+pw-2}" y="${M.t+13}" text-anchor="end"
        fill="${isZoomed ? '#a05c07' : '#9aa3b0'}"
        font-size="8.5" font-style="italic">${zoomTxt}</text>

  <!-- Interaction overlay (no duplicate id) -->
  <rect id="hover-overlay"
        x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="transparent" style="cursor:crosshair"/>

  <!-- Hover indicators -->
  <line id="h-line"
        x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t+ph}"
        stroke="#1a5fa8" stroke-width="1" stroke-dasharray="3,3"
        opacity="0" pointer-events="none"/>
  <circle id="h-dot" cx="${M.l}" cy="${M.t}" r="3.5"
          fill="#1a5fa8" stroke="white" stroke-width="1.2"
          opacity="0" pointer-events="none"/>
</svg>`;

  container.innerHTML = svg;
  _attachInteraction(container, curve, ord, geneNames,
                     { pw, ph, r0, r1, nVis, toX, toY });
}

// ── Interaction ───────────────────────────────────────────────
function _attachInteraction(container, curve, ord, geneNames,
                             { pw, ph, r0, r1, nVis, toX, toY }) {
  const svgEl   = container.querySelector('#es-svg');
  const overlay = container.querySelector('#hover-overlay');
  const hLine   = container.querySelector('#h-line');
  const hDot    = container.querySelector('#h-dot');
  const tooltip = document.getElementById('svg-tooltip');
  if (!svgEl || !overlay || !tooltip) return;

  // Cache CTM (screen px → SVG viewBox units).
  // Must be re-fetched if SVG is replaced (new _render call resets it).
  let _ctm = null;
  const getCTM   = () => { if (!_ctm) _ctm = svgEl.getScreenCTM(); return _ctm; };
  const toSVGPt  = (cx, cy) =>
    new DOMPoint(cx, cy).matrixTransform(getCTM().inverse());

  // Rank from SVG x coordinate (in viewBox units)
  const svgXtoRank = x =>
    Math.max(r0, Math.min(r1,
      r0 + Math.round((x - M.l) / pw * (nVis - 1))
    ));

  // ── Tooltip ────────────────────────────────────────────────
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

  // ── Scroll zoom ────────────────────────────────────────────
  overlay.addEventListener('wheel', e => {
    e.preventDefault();
    _ctm = null;   // invalidate after potential layout shift
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

  // ── Drag to pan ────────────────────────────────────────────
  // Key design: we do NOT re-render on every mousemove.
  // Instead we accumulate zoom state and re-render on mouseup.
  // This avoids the listener-removal bug caused by innerHTML replacement.
  //
  // For visual feedback during drag, we translate the SVG path
  // using a CSS transform (no DOM re-render required).

  let _drag = null;   // { startX (SVG), startI0, startI1, dx }

  overlay.addEventListener('mousedown', e => {
    if (e.button !== 0) return;
    e.preventDefault();
    _ctm = null;
    const pt = toSVGPt(e.clientX, e.clientY);
    _drag = {
      startSVGX: pt.x,
      startI0:   _zoom.i0,
      startI1:   _zoom.i1
    };
    _isDragging = false;
    overlay.style.cursor = 'grab';
  });

  // mousemove and mouseup are on window to handle fast mouse movement
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

    // dFrac > 0 means mouse moved right → show earlier genes → decrease i0
    const dFrac = -(dSVG / pw) * (_zoom.i1 - _zoom.i0 === 0
      ? 1 : (_drag.startI1 - _drag.startI0));
    // Recalculate relative to original span
    const span  = _drag.startI1 - _drag.startI0;
    const dFrac2 = -(dSVG / pw) * span;
    let i0 = _drag.startI0 + dFrac2;
    let i1 = _drag.startI1 + dFrac2;
    if (i0 < 0) { i0 = 0; i1 = span; }
    if (i1 > 1) { i1 = 1; i0 = 1 - span; }

    // Update zoom state without re-rendering (avoids listener removal)
    _zoom = { i0, i1 };

    // Visual feedback: translate the curve path using CSS transform
    // This is an approximation — full re-render on mouseup
    const esPath = container.querySelector('#es-path');
    if (esPath) {
      // dSVG in SVG units → pixel shift proportional to zoom
      const pxShift = dSVG * (pw / (span > 0 ? 1 : 1));
      // Keep it simple: just show a grabbing cursor, re-render on up
    }
  };

  const _onUp = e => {
    if (!_drag) return;
    const wasDragging = _isDragging;
    _drag        = null;
    _isDragging  = false;
    overlay.style.cursor = 'crosshair';

    if (wasDragging) {
      // Now re-render with updated _zoom
      // This replaces the SVG and re-attaches listeners correctly
      _render();
    }
  };

  window.addEventListener('mousemove', _onMove, { passive: true });
  window.addEventListener('mouseup',   _onUp);

  // Cleanup: when container content is replaced (next _render call),
  // remove the window-level listeners.
  const obs = new MutationObserver(() => {
    window.removeEventListener('mousemove', _onMove);
    window.removeEventListener('mouseup',   _onUp);
    obs.disconnect();
  });
  obs.observe(container, { childList: true });

  // ── Double-click reset ─────────────────────────────────────
  overlay.addEventListener('dblclick', () => {
    _zoom = { i0: 0, i1: 1 };
    _ctm  = null;
    _render();
  });

  // ── Resize: invalidate CTM ─────────────────────────────────
  // (listener is on window and will be cleaned up by MutationObserver)
  const _onResize = () => { _ctm = null; };
  window.addEventListener('resize', _onResize, { passive: true });
  // Register resize cleanup too
  const obs2 = new MutationObserver(() => {
    window.removeEventListener('resize', _onResize);
    obs2.disconnect();
  });
  obs2.observe(container, { childList: true });
}

// ── ES stats row ──────────────────────────────────────────────
export function renderESStats(result, engine, showFDR) {
  const fmt   = p => p == null ? '—'
    : Math.abs(p) < 0.001 ? p.toExponential(2) : p.toFixed(4);
  const sc    = p => p != null && p < 0.05 ? 'sig' : '';
  const isPar = engine === 'parametric';

  const cells = [
    { l:'ES',      v:result.es.toFixed(4),    c:result.es>=0?'pos':'neg' },
    { l:'NES',     v:result.nes.toFixed(3),   c:result.nes>=0?'pos':'neg' },
    { l:'NES-AD',  v:result.nes_ad.toFixed(3),c:'' },
    { l: isPar ? 'p<sub>KS</sub>&thinsp;(Γ)'  : 'p<sub>KS</sub>',
      v: fmt(result.pKS),     c: sc(result.pKS) },
    { l: isPar ? 'p<sub>AD</sub>&thinsp;(GΓ)' : 'p<sub>AD</sub>',
      v: fmt(result.pAD),     c: sc(result.pAD) },
    { l: 'p<sub>Cauchy</sub>',
      v: fmt(result.pCauchy), c: sc(result.pCauchy) },
    ...(showFDR && result.fdr != null
      ? [{ l:'FDR (BH)', v:fmt(result.fdr), c:sc(result.fdr) }] : []),
    { l:'Size', v:`${result.size}`, c:'' }
  ];

  document.getElementById('es-stats').innerHTML = cells.map(c =>
    `<div class="es-cell">
       <span class="es-lbl">${c.l}</span>
       <span class="es-val ${c.c}">${c.v}</span>
     </div>`
  ).join('');
}

const _f   = v => v.toFixed(2);
const _esc = s => (s??'').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
