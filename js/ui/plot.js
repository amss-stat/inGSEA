// ═══════════════════════════════════════════════════════════
//  ui/plot.js  ·  SVG enrichment curve with zoom
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

const M = { l: 64, r: 16, t: 22, b: 48, hitH: 24 };
const SVG_H        = 310;
const MAX_PATH_PTS = 4000;

// Zoom state per active result
let _zoom = { i0: 0, i1: 1 };   // normalised [0,1] over gene rank axis
let _currentResult = null;
let _currentPathways = null;
let _currentGenes = null;
let _currentContainer = null;

/**
 * Main entry: render enrichment curve SVG.
 */
export function drawCurve(result, pathways, geneNames, container) {
  _currentResult    = result;
  _currentPathways  = pathways;
  _currentGenes     = geneNames;
  _currentContainer = container;
  _renderWithZoom(result, pathways, geneNames, container, _zoom);
}

/** Reset zoom to full view. */
export function resetZoom() {
  _zoom = { i0: 0, i1: 1 };
  if (_currentResult)
    _renderWithZoom(_currentResult, _currentPathways,
                    _currentGenes, _currentContainer, _zoom);
}

// ── Core renderer (zoom-aware) ────────────────────────────────
function _renderWithZoom(result, pathways, geneNames, container, zoom) {
  const W  = Math.max(container.clientWidth || 680, 400);
  const pw = W - M.l - M.r;
  const ph = SVG_H - M.t - M.b - M.hitH - 6;

  const curve = result.curve;
  const nG    = curve.length;
  const ord   = result.obsOrd;
  const isPos = result.es >= 0;
  const lineC = isPos ? '#1c6e41' : '#b01c1c';
  const fillC = isPos ? 'rgba(28,110,65,.11)' : 'rgba(176,28,28,.09)';

  // Zoom window in gene-rank indices
  const r0 = Math.max(0,    Math.floor(zoom.i0 * (nG - 1)));
  const r1 = Math.min(nG-1, Math.ceil (zoom.i1 * (nG - 1)));
  const nVis = r1 - r0 + 1;  // number of visible genes

  // Y range over visible window
  let lo = 0, hi = 0;
  for (let i = r0; i <= r1; i++) {
    if (curve[i] < lo) lo = curve[i];
    if (curve[i] > hi) hi = curve[i];
  }
  const pad  = Math.max(-lo, hi) * 0.15 + 0.02;
  const yMin = lo - pad, yMax = hi + pad, yR = yMax - yMin;

  // Coordinate transforms (rank → SVG x, ES → SVG y)
  // toX maps visible rank r (r0..r1) to pixel x
  const toX = r => M.l + ((r - r0) / (nVis - 1)) * pw;
  const toY = es => M.t + ph * (1 - (es - yMin) / yR);
  const y0  = toY(0);

  // Subsampled path (always include peak if visible)
  const pkIdx = result.peakIdx;
  const step  = Math.max(1, Math.ceil(nVis / MAX_PATH_PTS));
  let pathD   = `M${_f(toX(r0))},${_f(toY(curve[r0]))}`;
  let fillD   = `M${_f(toX(r0))},${_f(y0)}`;

  for (let i = r0 + step; i <= r1; i += step) {
    // Always include peak point
    if (i > pkIdx && i - step <= pkIdx && pkIdx >= r0 && pkIdx <= r1) {
      pathD += ` L${_f(toX(pkIdx))},${_f(toY(curve[pkIdx]))}`;
      fillD += ` L${_f(toX(pkIdx))},${_f(toY(curve[pkIdx]))}`;
    }
    pathD += ` L${_f(toX(i))},${_f(toY(curve[i]))}`;
    fillD += ` L${_f(toX(i))},${_f(toY(curve[i]))}`;
  }
  // Ensure last point
  pathD += ` L${_f(toX(r1))},${_f(toY(curve[r1]))}`;
  fillD += ` L${_f(toX(r1))},${_f(toY(curve[r1]))} L${_f(toX(r1))},${_f(y0)} Z`;

  // Y ticks
  const nTicks = 5;
  const ticks  = Array.from({ length: nTicks + 1 },
    (_, k) => yMin + yR * k / nTicks);
  const tickSVG = ticks.map(v => {
    const y = _f(toY(v));
    return `<line x1="${M.l}" y1="${y}" x2="${M.l+pw}" y2="${y}"
      stroke="#e4e8ee" stroke-width="0.7" stroke-dasharray="3,4"/>
    <text x="${M.l-5}" y="${_f(toY(v)+3.5)}"
      text-anchor="end" fill="#8090a8" font-size="9.5">${v.toFixed(2)}</text>`;
  }).join('');

  // Hit strips (only visible window)
  const pathway = pathways.find(p => p.name === result.name);
  const mask    = pathway?.mask ?? null;
  const hitY    = M.t + ph + 6;
  let hitSVG    = '';
  if (mask) {
    const parts = [];
    for (let i = r0; i <= r1; i++) {
      if (mask[ord[i]]) parts.push(`M${_f(toX(i))},${hitY} v${M.hitH-5}`);
    }
    if (parts.length)
      hitSVG = `<path d="${parts.join(' ')}"
        stroke="${lineC}" stroke-width="1.1" opacity="0.65" fill="none"/>`;
  }

  // Peak marker (only if in visible window)
  let peakSVG = '';
  if (pkIdx >= r0 && pkIdx <= r1) {
    const pkX = _f(toX(pkIdx));
    const pkY = _f(toY(curve[pkIdx]));
    peakSVG = `
    <line x1="${pkX}" y1="${M.t}" x2="${pkX}" y2="${M.t+ph}"
          stroke="#a05c07" stroke-width="1" stroke-dasharray="4,3" opacity="0.7"
          clip-path="url(#cp)"/>
    <circle cx="${pkX}" cy="${pkY}" r="4"
            fill="#a05c07" stroke="white" stroke-width="1.2"
            clip-path="url(#cp)"/>`;
  }

  // X-axis labels (rank numbers, not normalised)
  const xMid   = _f(M.l + pw / 2);
  const xLabel = (r, anchor='middle') =>
    `<text x="${_f(toX(r))}" y="${hitY+M.hitH+14}"
      text-anchor="${anchor}" fill="#7a8698" font-size="9.5">${r+1}</text>`;

  // Zoom indicator
  const isZoomed = zoom.i0 > 0.001 || zoom.i1 < 0.999;
  const zoomBadge = isZoomed
    ? `<text x="${M.l+pw-2}" y="${M.t+13}"
         text-anchor="end" fill="#a05c07" font-size="9" font-weight="600">
         🔍 Rank ${r0+1}–${r1+1} · scroll to zoom · dbl-click to reset
       </text>`
    : `<text x="${M.l+pw-2}" y="${M.t+13}"
         text-anchor="end" fill="#9aa3b0" font-size="9">
         Scroll to zoom
       </text>`;

  const svg = `<svg xmlns="http://www.w3.org/2000/svg"
  viewBox="0 0 ${W} ${SVG_H}" width="${W}" height="${SVG_H}"
  style="font-family:'Inter',sans-serif;overflow:visible"
  id="es-svg" role="img" aria-label="Enrichment plot: ${_esc(result.name)}">
  <defs>
    <clipPath id="cp">
      <rect x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"/>
    </clipPath>
  </defs>
  <rect width="${W}" height="${SVG_H}" fill="#fafbfc"/>
  <rect x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="white" stroke="#dde2ea" stroke-width="0.8"/>
  ${tickSVG}
  ${[1,2,3,4].map(k=>{const x=_f(M.l+pw*k/5);
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
        fill="#9aa3b0" font-size="8" font-weight="600" letter-spacing=".06em">
    GENE HITS
  </text>
  <line x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t+ph}"
        stroke="#adb8c4" stroke-width="1.3"/>
  <line x1="${M.l}" y1="${M.t+ph}" x2="${M.l+pw}" y2="${M.t+ph}"
        stroke="#adb8c4" stroke-width="1.3"/>
  <text transform="translate(13,${_f(M.t+ph/2)}) rotate(-90)"
        text-anchor="middle" fill="#7a8698" font-size="10.5" font-weight="500">
    Enrichment Score
  </text>
  ${xLabel(r0, 'start')}
  <text x="${xMid}" y="${hitY+M.hitH+14}" text-anchor="middle"
        fill="#7a8698" font-size="9.5">${Math.floor((r0+r1)/2)+1}</text>
  ${xLabel(r1, 'end')}
  <text x="${xMid}" y="${hitY+M.hitH+30}" text-anchor="middle"
        fill="#7a8698" font-size="11" font-weight="500">Gene rank</text>
  ${zoomBadge}
  <!-- Hover + zoom overlay -->
  <rect id="hover-overlay" x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="transparent" style="cursor:crosshair"/>
  <line id="h-line" x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t+ph}"
        stroke="#1a5fa8" stroke-width="1" stroke-dasharray="3,3"
        opacity="0" pointer-events="none"/>
  <circle id="h-dot" cx="${M.l}" cy="${M.t}" r="3.5"
          fill="#1a5fa8" stroke="white" stroke-width="1.2"
          opacity="0" pointer-events="none"/>
</svg>`;

  container.innerHTML = svg;

  _attachInteraction(container, curve, ord, geneNames, nG,
                     r0, r1, nVis, toX, toY, result, pathways);
}

// ── Interaction: hover tooltip + scroll zoom ──────────────────
function _attachInteraction(container, curve, ord, geneNames, nG,
                             r0, r1, nVis, toX, toY,
                             result, pathways) {
  const overlay = container.querySelector('#hover-overlay');
  const hLine   = container.querySelector('#h-line');
  const hDot    = container.querySelector('#h-dot');
  const tooltip = document.getElementById('svg-tooltip');
  const svgEl   = container.querySelector('svg');
  if (!overlay || !tooltip) return;

  let ctm = null;

  // ── Tooltip (hover) ────────────────────────────────────────
  const getSVGPoint = (cx, cy) => {
    if (!ctm) ctm = svgEl.getScreenCTM();
    return new DOMPoint(cx, cy).matrixTransform(ctm.inverse());
  };

  overlay.addEventListener('mousemove', e => {
    const pt   = getSVGPoint(e.clientX, e.clientY);
    const frac = Math.max(0, Math.min(1, (pt.x - (container.getBoundingClientRect().left +
      /* M.l in SVG coords */ 0)) / 1));

    // Convert SVG x → rank within visible window
    const pw = svgEl.viewBox.baseVal.width - 64 - 16;  // M.l + M.r
    const svgX = pt.x;
    const rank = Math.max(r0, Math.min(r1,
      r0 + Math.round((svgX - 64) / pw * (nVis - 1))
    ));

    const gIdx = ord[rank];
    const name = geneNames[gIdx];
    const es   = curve[rank];
    const lx   = _f(toX(rank));
    const ly   = _f(toY(es));

    hLine.setAttribute('x1', lx); hLine.setAttribute('x2', lx);
    hLine.setAttribute('opacity', '0.55');
    hDot.setAttribute('cx', lx); hDot.setAttribute('cy', ly);
    hDot.setAttribute('opacity', '1');

    tooltip.style.display = 'block';
    tooltip.style.left    = (e.clientX + 14) + 'px';
    tooltip.style.top     = (e.clientY - 42) + 'px';
    tooltip.innerHTML     =
      `<strong>${_esc(name)}</strong><br>Rank #${rank+1}<br>ES&nbsp;${es.toFixed(4)}`;
  });

  overlay.addEventListener('mouseleave', () => {
    hLine.setAttribute('opacity', '0');
    hDot.setAttribute('opacity', '0');
    tooltip.style.display = 'none';
    ctm = null;
  });

  // ── Scroll zoom ────────────────────────────────────────────
  overlay.addEventListener('wheel', e => {
    e.preventDefault();

    // Mouse position as fraction of visible range [0,1]
    if (!ctm) ctm = svgEl.getScreenCTM();
    const pt  = getSVGPoint(e.clientX, e.clientY);
    const vb  = svgEl.viewBox.baseVal;
    const pw  = vb.width - 64 - 16;
    const foc = Math.max(0, Math.min(1, (pt.x - 64) / pw));

    const zoomFactor = e.deltaY < 0 ? 0.7 : 1 / 0.7;
    const curSpan    = _zoom.i1 - _zoom.i0;
    const newSpan    = Math.min(1, Math.max(0.01, curSpan * zoomFactor));

    // Keep focus point fixed
    let newI0 = _zoom.i0 + foc * (curSpan - newSpan);
    let newI1 = newI0 + newSpan;

    // Clamp to [0, 1]
    if (newI0 < 0) { newI0 = 0; newI1 = newSpan; }
    if (newI1 > 1) { newI1 = 1; newI0 = 1 - newSpan; }

    _zoom = { i0: newI0, i1: newI1 };
    ctm   = null;
    _renderWithZoom(result, pathways, geneNames, container, _zoom);
  }, { passive: false });

  // ── Double-click to reset zoom ─────────────────────────────
  overlay.addEventListener('dblclick', () => {
    _zoom = { i0: 0, i1: 1 };
    ctm   = null;
    _renderWithZoom(result, pathways, geneNames, container, _zoom);
  });

  window.addEventListener('resize', () => { ctm = null; }, { passive: true });
}

// ── Plot header & stats (unchanged API) ──────────────────────

export function updatePlotHeader(result) {
  document.getElementById('plot-name').textContent = result.name;
  document.getElementById('plot-section').style.display = 'block';
  const link = document.getElementById('plot-db-link');
  link.href  = msigdbUrl(result.name, result.url);
  link.style.display = 'inline-flex';
  // Reset zoom when switching pathways
  _zoom = { i0: 0, i1: 1 };
}

export function renderESStats(result, engine, showFDR) {
  const fmt = p => p == null ? '—'
    : Math.abs(p) < 0.001 ? p.toExponential(2) : p.toFixed(4);
  const sc = p => p != null && p < 0.05 ? 'sig' : '';
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
      ? [{ l:'FDR', v:fmt(result.fdr), c:sc(result.fdr) }] : []),
    { l:'Size', v:`${result.size}`, c:'' }
  ];

  document.getElementById('es-stats').innerHTML = cells.map(c => `
    <div class="es-cell">
      <span class="es-lbl">${c.l}</span>
      <span class="es-val ${c.c}">${c.v}</span>
    </div>`).join('');
}

const _f   = v => v.toFixed(2);
const _esc = s => (s??'').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
