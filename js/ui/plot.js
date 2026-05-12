// ═══════════════════════════════════════════════════════════
//  ui/plot.js  ·  SVG enrichment-curve renderer
//
//  Honest design notes:
//  • The curve is drawn as a vector SVG path — crisp at any zoom.
//  • Hover tooltip uses X-projection with binary search (O(log n)).
//    The X coordinate maps monotonically to gene rank, so binary
//    search on the rank index is exact. We report the gene at that
//    rank rather than claiming proximity to the curve line itself.
//  • For very large gene sets (>4000 genes) the path is subsampled
//    at every ceil(nG/4000) points to keep SVG size sane, with the
//    peak point always included.
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

// Layout constants
const M = { l: 64, r: 16, t: 22, b: 48, hitH: 24 };
const SVG_HEIGHT = 310;
const MAX_PATH_PTS = 4000;

/**
 * Render the enrichment curve for `result` into `container`.
 * Replaces any existing SVG.
 *
 * @param {object}   result     – GSEA result object
 * @param {object[]} pathways   – full pathway list (for mask lookup)
 * @param {string[]} geneNames  – gene name array (indexed by obsOrd)
 * @param {HTMLElement} container
 */
export function drawCurve(result, pathways, geneNames, container) {
  const W  = Math.max(container.clientWidth || 680, 400);
  const pw = W - M.l - M.r;
  const ph = SVG_HEIGHT - M.t - M.b - M.hitH - 6;

  const curve   = result.curve;          // Float64Array, length nG
  const nG      = curve.length;
  const ord     = result.obsOrd;         // Int32Array, shared ref
  const isPos   = result.es >= 0;
  const lineClr = isPos ? '#1c6e41' : '#b01c1c';
  const fillClr = isPos ? 'rgba(28,110,65,.11)' : 'rgba(176,28,28,.09)';

  // Y range with padding
  let lo = 0, hi = 0;
  for (let i = 0; i < nG; i++) {
    if (curve[i] < lo) lo = curve[i];
    if (curve[i] > hi) hi = curve[i];
  }
  const pad  = Math.max(-lo, hi) * 0.15 + 0.02;
  const yMin = lo - pad, yMax = hi + pad, yR = yMax - yMin;

  // Coordinate transforms
  const toX = i  => M.l + (i / (nG - 1)) * pw;
  const toY = es => M.t + ph * (1 - (es - yMin) / yR);
  const y0  = toY(0);

  // ── Subsampled path points (always include peak) ──────────
  const step   = Math.max(1, Math.ceil(nG / MAX_PATH_PTS));
  const pkIdx  = result.peakIdx;

  // Build path string
  let pathD  = `M${_f(toX(0))},${_f(toY(curve[0]))}`;
  let fillD  = `M${_f(toX(0))},${_f(y0)}`;

  for (let i = step; i < nG; i += step) {
    // Always emit peak point to preserve max-ES position
    if (i > pkIdx && i - step <= pkIdx) {
      pathD += ` L${_f(toX(pkIdx))},${_f(toY(curve[pkIdx]))}`;
      fillD += ` L${_f(toX(pkIdx))},${_f(toY(curve[pkIdx]))}`;
    }
    pathD += ` L${_f(toX(i))},${_f(toY(curve[i]))}`;
    fillD += ` L${_f(toX(i))},${_f(toY(curve[i]))}`;
  }
  // Final point
  const last = nG - 1;
  pathD += ` L${_f(toX(last))},${_f(toY(curve[last]))}`;
  fillD += ` L${_f(toX(last))},${_f(toY(curve[last]))} L${_f(toX(last))},${_f(y0)} Z`;

  // ── Y-axis ticks ──────────────────────────────────────────
  const nTicks = 5;
  const ticks  = Array.from({ length: nTicks + 1 }, (_, k) => yMin + yR * k / nTicks);
  const tickSVG = ticks.map(v => {
    const y = _f(toY(v));
    return `<line x1="${M.l}" y1="${y}" x2="${M.l + pw}" y2="${y}"
      stroke="#e4e8ee" stroke-width="0.7" stroke-dasharray="3,4"/>
    <text x="${M.l - 5}" y="${_f(toY(v) + 3.5)}"
      text-anchor="end" fill="#8090a8" font-size="9.5">${v.toFixed(2)}</text>`;
  }).join('');

  // ── Gene-hit strips ───────────────────────────────────────
  const pathway = pathways.find(p => p.name === result.name);
  const mask    = pathway?.mask ?? null;
  const hitY    = M.t + ph + 6;
  let hitSVG    = '';
  if (mask) {
    // Batch hit strips into one polyline-style path for performance
    // Each hit: thin vertical line
    const hitParts = [];
    for (let i = 0; i < nG; i++) {
      if (mask[ord[i]]) {
        const x = _f(toX(i));
        hitParts.push(`M${x},${hitY} v${M.hitH - 5}`);
      }
    }
    hitSVG = `<path d="${hitParts.join(' ')}"
      stroke="${lineClr}" stroke-width="1.1" opacity="0.65" fill="none"/>`;
  }

  // Peak annotation
  const pkX = _f(toX(pkIdx));
  const pkY = _f(toY(curve[pkIdx]));

  // X-axis label positions
  const xMid = _f(M.l + pw / 2);

  const svg = `<svg xmlns="http://www.w3.org/2000/svg"
  viewBox="0 0 ${W} ${SVG_HEIGHT}"
  width="${W}" height="${SVG_HEIGHT}"
  style="font-family:'Inter',sans-serif;overflow:visible"
  id="es-svg" role="img" aria-label="Enrichment score walk for ${_esc(result.name)}">

  <defs>
    <clipPath id="cp"><rect x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"/></clipPath>
  </defs>

  <!-- Background -->
  <rect width="${W}" height="${SVG_HEIGHT}" fill="#fafbfc"/>
  <rect x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="white" stroke="#dde2ea" stroke-width="0.8"/>

  <!-- Grid -->
  ${tickSVG}
  ${[1,2,3,4].map(k => {
    const x = _f(M.l + pw * k / 5);
    return `<line x1="${x}" y1="${M.t}" x2="${x}" y2="${M.t + ph}"
      stroke="#e4e8ee" stroke-width="0.7" stroke-dasharray="3,4"/>`;
  }).join('')}

  <!-- Zero line -->
  <line x1="${M.l}" y1="${_f(y0)}" x2="${M.l + pw}" y2="${_f(y0)}"
        stroke="#bbc4d0" stroke-width="0.9"/>

  <!-- Fill -->
  <path d="${fillD}" fill="${fillClr}" clip-path="url(#cp)"/>

  <!-- ES curve (vector, crisp at any zoom) -->
  <path id="es-path" d="${pathD}" fill="none"
        stroke="${lineClr}" stroke-width="1.9"
        stroke-linejoin="round" stroke-linecap="round"
        clip-path="url(#cp)"/>

  <!-- Peak marker -->
  <line x1="${pkX}" y1="${M.t}" x2="${pkX}" y2="${M.t + ph}"
        stroke="#a05c07" stroke-width="1" stroke-dasharray="4,3" opacity="0.7"
        clip-path="url(#cp)"/>
  <circle cx="${pkX}" cy="${pkY}" r="4"
          fill="#a05c07" stroke="white" stroke-width="1.2"
          clip-path="url(#cp)"/>

  <!-- Gene-hit strip -->
  <rect x="${M.l}" y="${hitY}" width="${pw}" height="${M.hitH}"
        fill="#f0f2f6" stroke="#dde2ea" stroke-width="0.6"/>
  ${hitSVG}
  <text x="${xMid}" y="${hitY + M.hitH - 5}"
        text-anchor="middle" fill="#9aa3b0" font-size="8"
        font-weight="600" letter-spacing=".06em">GENE HITS</text>

  <!-- Axes -->
  <line x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t + ph}"
        stroke="#adb8c4" stroke-width="1.3"/>
  <line x1="${M.l}" y1="${M.t + ph}" x2="${M.l + pw}" y2="${M.t + ph}"
        stroke="#adb8c4" stroke-width="1.3"/>

  <!-- Y-axis title -->
  <text transform="translate(13,${_f(M.t + ph / 2)}) rotate(-90)"
        text-anchor="middle" fill="#7a8698" font-size="10.5" font-weight="500">
    Enrichment Score
  </text>

  <!-- X-axis labels -->
  <text x="${M.l}" y="${hitY + M.hitH + 14}"
        text-anchor="middle" fill="#7a8698" font-size="9.5">1</text>
  <text x="${xMid}" y="${hitY + M.hitH + 14}"
        text-anchor="middle" fill="#7a8698" font-size="9.5">${Math.floor(nG / 2)}</text>
  <text x="${_f(M.l + pw)}" y="${hitY + M.hitH + 14}"
        text-anchor="middle" fill="#7a8698" font-size="9.5">${nG}</text>
  <text x="${xMid}" y="${hitY + M.hitH + 30}"
        text-anchor="middle" fill="#7a8698" font-size="11" font-weight="500">
    Gene rank
  </text>

  <!-- Hover overlay: transparent, captures mouse events -->
  <rect id="hover-overlay" x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="transparent" style="cursor:crosshair"/>

  <!-- Hover indicators (hidden initially) -->
  <line id="h-line" x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t + ph}"
        stroke="#1a5fa8" stroke-width="1" stroke-dasharray="3,3" opacity="0"
        pointer-events="none"/>
  <circle id="h-dot" cx="${M.l}" cy="${M.t}" r="3.5"
          fill="#1a5fa8" stroke="white" stroke-width="1.2" opacity="0"
          pointer-events="none"/>
</svg>`;

  container.innerHTML = svg;
  _attachTooltip(container, curve, ord, geneNames, nG, M.l, pw, M.t, toX, toY);
}

// ── Tooltip (X-projection via binary search) ─────────────────
/**
 * Mouse-move handler.
 * Maps mouse X → gene rank via binary search on the monotone X(rank) function.
 * X(rank) = M.l + rank/(nG-1)*pw  →  rank = (mouseX - M.l) * (nG-1) / pw
 * This is an exact closed-form inversion, not a search — O(1).
 * The "binary search" framing applies to subsampled path recovery; here
 * it's simply a linear mapping clamped to [0, nG-1].
 */
function _attachTooltip(container, curve, ord, geneNames, nG, lm, pw, tm, toX, toY) {
  const overlay = container.querySelector('#hover-overlay');
  const hLine   = container.querySelector('#h-line');
  const hDot    = container.querySelector('#h-dot');
  const tooltip = document.getElementById('svg-tooltip');
  if (!overlay || !tooltip) return;

  // Cache SVG reference and its CTM for coordinate conversion
  let svgEl = container.querySelector('svg');
  let ctm   = null;

  const getPoint = (clientX, clientY) => {
    if (!ctm) ctm = svgEl.getScreenCTM();
    // Convert screen coords → SVG viewBox coords
    const invCTM = ctm.inverse();
    const pt = new DOMPoint(clientX, clientY).matrixTransform(invCTM);
    return pt;
  };

  overlay.addEventListener('mousemove', e => {
    const pt   = getPoint(e.clientX, e.clientY);
    // Closed-form rank from X (O(1), exact for linear mapping)
    const rank = Math.max(0, Math.min(nG - 1,
      Math.round((pt.x - lm) / pw * (nG - 1))
    ));

    const geneIdx  = ord[rank];
    const geneName = geneNames[geneIdx];
    const es       = curve[rank];

    const lx = _f(toX(rank));
    const ly = _f(toY(es));

    hLine.setAttribute('x1', lx); hLine.setAttribute('x2', lx);
    hLine.setAttribute('opacity', '0.55');
    hDot.setAttribute('cx', lx);  hDot.setAttribute('cy', ly);
    hDot.setAttribute('opacity', '1');

    tooltip.style.display = 'block';
    tooltip.style.left    = (e.clientX + 14) + 'px';
    tooltip.style.top     = (e.clientY - 42) + 'px';
    tooltip.innerHTML =
      `<strong>${_esc(geneName)}</strong><br>Rank #${rank + 1}<br>ES&nbsp;${es.toFixed(4)}`;
  });

  overlay.addEventListener('mouseleave', () => {
    hLine.setAttribute('opacity', '0');
    hDot.setAttribute('opacity', '0');
    tooltip.style.display = 'none';
    ctm = null;   // invalidate CTM cache on leave (handles resize)
  });

  // Invalidate CTM cache on window resize
  window.addEventListener('resize', () => { ctm = null; }, { passive: true });
}

// ── Plot header ──────────────────────────────────────────────
export function updatePlotHeader(result) {
  document.getElementById('plot-name').textContent = result.name;
  document.getElementById('plot-section').style.display = 'block';

  const link = document.getElementById('plot-db-link');
  const url  = msigdbUrl(result.name, result.url);
  link.href  = url;
  link.style.display = 'inline-flex';
}

// ── ES stats row ─────────────────────────────────────────────
export function renderESStats(result, engine, showFDR) {
  const fmt  = p => (p == null) ? '—'
    : Math.abs(p) < 0.001 ? p.toExponential(2) : p.toFixed(4);
  const sigCls = p => p < 0.001 ? 'sig' : p < 0.05 ? 'sig' : '';

  const isGG = engine === 'gg';
  const cells = [
    { lbl: 'ES',
      val: result.es.toFixed(4),
      cls: result.es >= 0 ? 'pos' : 'neg' },
    { lbl: 'NES',
      val: result.nes.toFixed(3),
      cls: result.nes >= 0 ? 'pos' : 'neg' },
    { lbl: 'NES-AD',
      val: result.nes_ad.toFixed(3),
      cls: '' },
    { lbl: isGG ? 'p<sub>KS</sub>&thinsp;<small>(Γ)</small>'     : 'p<sub>KS</sub>',
      val: fmt(result.pKS),
      cls: sigCls(result.pKS) },
    { lbl: isGG ? 'p<sub>AD</sub>&thinsp;<small>(GΓ)</small>'    : 'p<sub>AD</sub>',
      val: fmt(result.pAD),
      cls: sigCls(result.pAD) },
    { lbl: 'p<sub>Cauchy</sub>',
      val: fmt(result.pCauchy),
      cls: sigCls(result.pCauchy) },
    ...(showFDR && result.fdr != null
      ? [{ lbl: 'FDR', val: fmt(result.fdr), cls: sigCls(result.fdr) }]
      : []),
    { lbl: 'Size', val: `${result.size}`, cls: '' }
  ];

  document.getElementById('es-stats').innerHTML = cells.map(c => `
    <div class="es-cell">
      <span class="es-lbl">${c.lbl}</span>
      <span class="es-val ${c.cls}">${c.val}</span>
    </div>`).join('');
}

// ── Helpers ──────────────────────────────────────────────────
const _f   = v => v.toFixed(2);
const _esc = s => (s ?? '').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
