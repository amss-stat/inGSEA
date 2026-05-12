// ═══════════════════════════════════════════════════════════
//  ui/plot.js  ·  v2.4
//
//  Changes:
//  • Fixed hover coordinate calculation (use SVG viewBox coords)
//  • Drag-to-pan after zoom (mouse drag on plot area)
//  • "Reset zoom" button rendered in SVG
//  • Zoom state resets on pathway switch
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

const M          = { l: 64, r: 16, t: 22, b: 48, hitH: 24 };
const SVG_H      = 310;
const MAX_PTS    = 4000;

// Module-level zoom state
let _zoom = { i0: 0, i1: 1 };
let _cur  = { result: null, pathways: null, genes: null, container: null };

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
  _zoom = { i0: 0, i1: 1 };   // reset zoom on pathway change
}

// ── Core renderer ─────────────────────────────────────────────
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

  // Visible gene rank window
  const r0   = Math.max(0,    Math.floor(_zoom.i0 * (nG - 1)));
  const r1   = Math.min(nG-1, Math.ceil (_zoom.i1 * (nG - 1)));
  const nVis = r1 - r0 + 1;

  // Y range over visible window
  let lo = 0, hi = 0;
  for (let i = r0; i <= r1; i++) {
    if (curve[i] < lo) lo = curve[i];
    if (curve[i] > hi) hi = curve[i];
  }
  const pad  = Math.max(-lo, hi) * 0.15 + 0.02;
  const yMin = lo - pad, yMax = hi + pad, yR = yMax - yMin;

  // Coordinate transforms (in SVG viewBox units, not pixels)
  // toX(rank): maps rank index r (r0..r1) → SVG x coordinate
  const toX = r  => M.l + ((r - r0) / Math.max(nVis - 1, 1)) * pw;
  const toY = es => M.t + ph * (1 - (es - yMin) / yR);
  const y0  = toY(0);

  // Path with subsampling
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

  // Y ticks
  const ticks = Array.from({ length: 6 }, (_, k) => yMin + yR * k / 5);
  const tickSVG = ticks.map(v => {
    const y = _f(toY(v));
    return `<line x1="${M.l}" y1="${y}" x2="${M.l+pw}" y2="${y}"
      stroke="#e4e8ee" stroke-width="0.7" stroke-dasharray="3,4"/>
    <text x="${M.l-5}" y="${_f(toY(v)+3.5)}"
      text-anchor="end" fill="#8090a8" font-size="9.5">${v.toFixed(2)}</text>`;
  }).join('');

  // Hit strips (visible window)
  const pathway = pathways.find(p => p.name === result.name);
  const mask    = pathway?.mask ?? null;
  const hitY    = M.t + ph + 6;
  let hitSVG    = '';
  if (mask) {
    const parts = [];
    for (let i = r0; i <= r1; i++)
      if (mask[ord[i]]) parts.push(`M${_f(toX(i))},${hitY} v${M.hitH-5}`);
    if (parts.length)
      hitSVG = `<path d="${parts.join(' ')}"
        stroke="${lineC}" stroke-width="1.1" opacity="0.65" fill="none"/>`;
  }

  // Peak (if visible)
  let peakSVG = '';
  if (pkIdx >= r0 && pkIdx <= r1) {
    const px = _f(toX(pkIdx)), py = _f(toY(curve[pkIdx]));
    peakSVG = `<line x1="${px}" y1="${M.t}" x2="${px}" y2="${M.t+ph}"
      stroke="#a05c07" stroke-width="1" stroke-dasharray="4,3" opacity="0.7"
      clip-path="url(#cp)"/>
    <circle cx="${px}" cy="${py}" r="4"
      fill="#a05c07" stroke="white" stroke-width="1.2" clip-path="url(#cp)"/>`;
  }

  const xMid     = _f(M.l + pw / 2);
  const isZoomed = _zoom.i0 > 0.001 || _zoom.i1 < 0.999;
  const zoomInfo = isZoomed
    ? `Rank ${r0+1}–${r1+1} · Drag to pan · Scroll to zoom · Dbl-click reset`
    : `Scroll to zoom · Drag to pan`;

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
    const x = _f(M.l + pw*k/5);
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
  <text x="${_f(toX(r0))}" y="${hitY+M.hitH+14}"
        text-anchor="start" fill="#7a8698" font-size="9.5">${r0+1}</text>
  <text x="${xMid}" y="${hitY+M.hitH+14}"
        text-anchor="middle" fill="#7a8698" font-size="9.5">
    ${Math.round((r0+r1)/2)+1}
  </text>
  <text x="${_f(toX(r1))}" y="${hitY+M.hitH+14}"
        text-anchor="end" fill="#7a8698" font-size="9.5">${r1+1}</text>
  <text x="${xMid}" y="${hitY+M.hitH+30}"
        text-anchor="middle" fill="#7a8698" font-size="11" font-weight="500">
    Gene rank
  </text>
  <text x="${M.l+pw-2}" y="${M.t+13}"
        text-anchor="end" fill="${isZoomed?'#a05c07':'#9aa3b0'}"
        font-size="8.5" font-style="italic">${zoomInfo}</text>
  <!-- Interaction overlay -->
  <rect id="hover-overlay" x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="transparent" style="cursor:crosshair" id="hover-overlay"/>
  <line id="h-line" x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t+ph}"
        stroke="#1a5fa8" stroke-width="1" stroke-dasharray="3,3"
        opacity="0" pointer-events="none"/>
  <circle id="h-dot" cx="${M.l}" cy="${M.t}" r="3.5"
          fill="#1a5fa8" stroke="white" stroke-width="1.2"
          opacity="0" pointer-events="none"/>
</svg>`;

  container.innerHTML = svg;

  // Store layout for event handlers
  const layout = { W, pw, ph, r0, r1, nVis, nG, toX, toY, M };
  _attachInteraction(container, curve, ord, geneNames, layout, result, pathways);
}

// ── Interaction handler ───────────────────────────────────────
function _attachInteraction(container, curve, ord, geneNames, layout, result, pathways) {
  const overlay = container.querySelector('#hover-overlay');
  const hLine   = container.querySelector('#h-line');
  const hDot    = container.querySelector('#h-dot');
  const tooltip = document.getElementById('svg-tooltip');
  const svgEl   = container.querySelector('#es-svg');
  if (!overlay || !tooltip || !svgEl) return;

  const { pw, ph, r0, r1, nVis, nG, toX, toY, M: m } = layout;

  // CTM: screen → SVG viewBox coordinate transform
  // We cache and invalidate on resize / scroll
  let _ctm = null;
  const getCTM = () => {
    if (!_ctm) _ctm = svgEl.getScreenCTM();
    return _ctm;
  };
  const toSVG = (cx, cy) => {
    const pt = new DOMPoint(cx, cy).matrixTransform(getCTM().inverse());
    return pt;
  };

  // Convert SVG x → gene rank (in visible window)
  // All coordinates here are in SVG viewBox units (not pixels)
  const svgXtoRank = svgX =>
    Math.max(r0, Math.min(r1,
      r0 + Math.round((svgX - m.l) / pw * (nVis - 1))
    ));

  // ── Tooltip ────────────────────────────────────────────────
  let _dragging = false;

  overlay.addEventListener('mousemove', e => {
    if (_dragging) return;
    const pt   = toSVG(e.clientX, e.clientY);
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
    if (_dragging) return;
    hLine.setAttribute('opacity', '0');
    hDot.setAttribute('opacity', '0');
    tooltip.style.display = 'none';
    _ctm = null;
  });

  // ── Scroll zoom ────────────────────────────────────────────
  overlay.addEventListener('wheel', e => {
    e.preventDefault();
    const pt       = toSVG(e.clientX, e.clientY);
    const foc      = Math.max(0, Math.min(1, (pt.x - m.l) / pw));
    const zFactor  = e.deltaY < 0 ? 0.7 : 1 / 0.7;
    const curSpan  = _zoom.i1 - _zoom.i0;
    const newSpan  = Math.max(0.01, Math.min(1, curSpan * zFactor));
    let   newI0    = _zoom.i0 + foc * (curSpan - newSpan);
    let   newI1    = newI0 + newSpan;
    if (newI0 < 0) { newI0 = 0; newI1 = newSpan; }
    if (newI1 > 1) { newI1 = 1; newI0 = 1 - newSpan; }
    _zoom = { i0: newI0, i1: newI1 };
    _ctm  = null;
    _render();
  }, { passive: false });

  // ── Drag to pan ────────────────────────────────────────────
  let _dragStart = null;   // { svgX, i0, i1 }

  overlay.addEventListener('mousedown', e => {
    if (e.button !== 0) return;
    const pt    = toSVG(e.clientX, e.clientY);
    _dragStart  = { svgX: pt.x, i0: _zoom.i0, i1: _zoom.i1 };
    _dragging   = false;
    overlay.style.cursor = 'grabbing';
    e.preventDefault();
  });

  // Attach move/up to window so drag works outside SVG
  const onMove = e => {
    if (!_dragStart) return;
    const pt    = toSVG(e.clientX, e.clientY);
    const dSVG  = pt.x - _dragStart.svgX;   // delta in SVG units
    const dFrac = -dSVG / pw;                // fraction of full range
    const span  = _dragStart.i1 - _dragStart.i0;
    let newI0   = _dragStart.i0 + dFrac;
    let newI1   = _dragStart.i1 + dFrac;
    if (newI0 < 0) { newI0 = 0; newI1 = span; }
    if (newI1 > 1) { newI1 = 1; newI0 = 1 - span; }

    // Only trigger a re-render if we've actually moved
    if (Math.abs(dSVG) > 2) {
      _dragging = true;
      hLine.setAttribute('opacity', '0');
      hDot.setAttribute('opacity', '0');
      tooltip.style.display = 'none';
    }
    if (_dragging) {
      _zoom = { i0: newI0, i1: newI1 };
      _ctm  = null;
      _render();
    }
  };

  const onUp = () => {
    if (_dragStart) {
      _dragStart  = null;
      _dragging   = false;
      overlay.style.cursor = 'crosshair';
    }
  };

  window.addEventListener('mousemove', onMove);
  window.addEventListener('mouseup',   onUp);

  // Cleanup on next drawCurve call (container.innerHTML replacement removes SVG)
  const observer = new MutationObserver(() => {
    window.removeEventListener('mousemove', onMove);
    window.removeEventListener('mouseup',   onUp);
    observer.disconnect();
  });
  observer.observe(container, { childList: true });

  // ── Double-click reset ─────────────────────────────────────
  overlay.addEventListener('dblclick', () => {
    _zoom = { i0: 0, i1: 1 };
    _ctm  = null;
    _render();
  });

  window.addEventListener('resize', () => { _ctm = null; }, { passive: true });
}

// ── ES stats row ──────────────────────────────────────────────
export function renderESStats(result, engine, showFDR) {
  const fmt = p => p == null ? '—'
    : Math.abs(p) < 0.001 ? p.toExponential(2) : p.toFixed(4);
  const sc  = p => p != null && p < 0.05 ? 'sig' : '';
  const isPar = engine === 'parametric';

  const cells = [
    { l: 'ES',      v: result.es.toFixed(4),     c: result.es>=0?'pos':'neg' },
    { l: 'NES',     v: result.nes.toFixed(3),    c: result.nes>=0?'pos':'neg' },
    { l: 'NES-AD',  v: result.nes_ad.toFixed(3), c: '' },
    { l: isPar ? 'p<sub>KS</sub>&thinsp;(Γ)'  : 'p<sub>KS</sub>',
      v: fmt(result.pKS),     c: sc(result.pKS) },
    { l: isPar ? 'p<sub>AD</sub>&thinsp;(GΓ)' : 'p<sub>AD</sub>',
      v: fmt(result.pAD),     c: sc(result.pAD) },
    { l: 'p<sub>Cauchy</sub>',
      v: fmt(result.pCauchy), c: sc(result.pCauchy) },
    ...(showFDR && result.fdr != null
      ? [{ l:'FDR (BH)', v:fmt(result.fdr), c:sc(result.fdr) }] : []),
    { l: 'Size',    v: `${result.size}`,          c: '' }
  ];

  document.getElementById('es-stats').innerHTML = cells.map(c => `
    <div class="es-cell">
      <span class="es-lbl">${c.l}</span>
      <span class="es-val ${c.c}">${c.v}</span>
    </div>`).join('');
}

const _f   = v => v.toFixed(2);
const _esc = s => (s ?? '').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
