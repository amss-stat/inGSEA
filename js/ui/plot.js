// ═══════════════════════════════════════════════════════════
//  ui/plot.js  ·  SVG enrichment curve
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

const M = { l: 64, r: 16, t: 22, b: 48, hitH: 24 };
const SVG_HEIGHT     = 310;
const MAX_PATH_PTS   = 4000;

/**
 * Render enrichment curve as inline SVG.
 */
export function drawCurve(result, pathways, geneNames, container) {
  const W  = Math.max(container.clientWidth || 680, 400);
  const pw = W - M.l - M.r;
  const ph = SVG_HEIGHT - M.t - M.b - M.hitH - 6;

  const curve = result.curve;
  const nG    = curve.length;
  const ord   = result.obsOrd;
  const isPos = result.es >= 0;
  const lineC = isPos ? '#1c6e41' : '#b01c1c';
  const fillC = isPos ? 'rgba(28,110,65,.11)' : 'rgba(176,28,28,.09)';

  let lo = 0, hi = 0;
  for (let i = 0; i < nG; i++) {
    if (curve[i] < lo) lo = curve[i];
    if (curve[i] > hi) hi = curve[i];
  }
  const pad  = Math.max(-lo, hi) * 0.15 + 0.02;
  const yMin = lo - pad, yMax = hi + pad, yR = yMax - yMin;

  const toX = i  => M.l + (i / (nG - 1)) * pw;
  const toY = es => M.t + ph * (1 - (es - yMin) / yR);
  const y0  = toY(0);

  // Subsampled path
  const step  = Math.max(1, Math.ceil(nG / MAX_PATH_PTS));
  const pkIdx = result.peakIdx;

  let pathD = `M${_f(toX(0))},${_f(toY(curve[0]))}`;
  let fillD = `M${_f(toX(0))},${_f(y0)}`;

  for (let i = step; i < nG; i += step) {
    if (i > pkIdx && i - step <= pkIdx) {
      pathD += ` L${_f(toX(pkIdx))},${_f(toY(curve[pkIdx]))}`;
      fillD += ` L${_f(toX(pkIdx))},${_f(toY(curve[pkIdx]))}`;
    }
    pathD += ` L${_f(toX(i))},${_f(toY(curve[i]))}`;
    fillD += ` L${_f(toX(i))},${_f(toY(curve[i]))}`;
  }
  const last = nG - 1;
  pathD += ` L${_f(toX(last))},${_f(toY(curve[last]))}`;
  fillD += ` L${_f(toX(last))},${_f(toY(curve[last]))} L${_f(toX(last))},${_f(y0)} Z`;

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

  // Hit strips
  const pathway = pathways.find(p => p.name === result.name);
  const mask    = pathway?.mask ?? null;
  const hitY    = M.t + ph + 6;
  let hitSVG    = '';
  if (mask) {
    const parts = [];
    for (let i = 0; i < nG; i++) {
      if (mask[ord[i]]) parts.push(`M${_f(toX(i))},${hitY} v${M.hitH-5}`);
    }
    hitSVG = `<path d="${parts.join(' ')}"
      stroke="${lineC}" stroke-width="1.1" opacity="0.65" fill="none"/>`;
  }

  const pkX  = _f(toX(pkIdx));
  const pkY  = _f(toY(curve[pkIdx]));
  const xMid = _f(M.l + pw / 2);

  const svg = `<svg xmlns="http://www.w3.org/2000/svg"
  viewBox="0 0 ${W} ${SVG_HEIGHT}" width="${W}" height="${SVG_HEIGHT}"
  style="font-family:'Inter',sans-serif;overflow:visible"
  id="es-svg" role="img" aria-label="Enrichment plot: ${_esc(result.name)}">
  <defs>
    <clipPath id="cp"><rect x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"/></clipPath>
  </defs>
  <rect width="${W}" height="${SVG_HEIGHT}" fill="#fafbfc"/>
  <rect x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="white" stroke="#dde2ea" stroke-width="0.8"/>
  ${tickSVG}
  ${[1,2,3,4].map(k=>{const x=_f(M.l+pw*k/5);return`<line x1="${x}" y1="${M.t}" x2="${x}" y2="${M.t+ph}" stroke="#e4e8ee" stroke-width="0.7" stroke-dasharray="3,4"/>`;}).join('')}
  <line x1="${M.l}" y1="${_f(y0)}" x2="${M.l+pw}" y2="${_f(y0)}"
        stroke="#bbc4d0" stroke-width="0.9"/>
  <path d="${fillD}" fill="${fillC}" clip-path="url(#cp)"/>
  <path id="es-path" d="${pathD}" fill="none" stroke="${lineC}" stroke-width="1.9"
        stroke-linejoin="round" stroke-linecap="round" clip-path="url(#cp)"/>
  <line x1="${pkX}" y1="${M.t}" x2="${pkX}" y2="${M.t+ph}"
        stroke="#a05c07" stroke-width="1" stroke-dasharray="4,3" opacity="0.7" clip-path="url(#cp)"/>
  <circle cx="${pkX}" cy="${pkY}" r="4" fill="#a05c07" stroke="white" stroke-width="1.2" clip-path="url(#cp)"/>
  <rect x="${M.l}" y="${hitY}" width="${pw}" height="${M.hitH}"
        fill="#f0f2f6" stroke="#dde2ea" stroke-width="0.6"/>
  ${hitSVG}
  <text x="${xMid}" y="${hitY+M.hitH-5}" text-anchor="middle" fill="#9aa3b0"
        font-size="8" font-weight="600" letter-spacing=".06em">GENE HITS</text>
  <line x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t+ph}" stroke="#adb8c4" stroke-width="1.3"/>
  <line x1="${M.l}" y1="${M.t+ph}" x2="${M.l+pw}" y2="${M.t+ph}" stroke="#adb8c4" stroke-width="1.3"/>
  <text transform="translate(13,${_f(M.t+ph/2)}) rotate(-90)"
        text-anchor="middle" fill="#7a8698" font-size="10.5" font-weight="500">Enrichment Score</text>
  <text x="${M.l}" y="${hitY+M.hitH+14}" text-anchor="middle" fill="#7a8698" font-size="9.5">1</text>
  <text x="${xMid}" y="${hitY+M.hitH+14}" text-anchor="middle" fill="#7a8698" font-size="9.5">${Math.floor(nG/2)}</text>
  <text x="${_f(M.l+pw)}" y="${hitY+M.hitH+14}" text-anchor="middle" fill="#7a8698" font-size="9.5">${nG}</text>
  <text x="${xMid}" y="${hitY+M.hitH+30}" text-anchor="middle" fill="#7a8698" font-size="11" font-weight="500">Gene rank</text>
  <rect id="hover-overlay" x="${M.l}" y="${M.t}" width="${pw}" height="${ph}"
        fill="transparent" style="cursor:crosshair"/>
  <line id="h-line" x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${M.t+ph}"
        stroke="#1a5fa8" stroke-width="1" stroke-dasharray="3,3" opacity="0" pointer-events="none"/>
  <circle id="h-dot" cx="${M.l}" cy="${M.t}" r="3.5"
          fill="#1a5fa8" stroke="white" stroke-width="1.2" opacity="0" pointer-events="none"/>
</svg>`;

  container.innerHTML = svg;
  _attachTooltip(container, curve, ord, geneNames, nG, M.l, pw, toX, toY);
}

function _attachTooltip(container, curve, ord, geneNames, nG, lm, pw, toX, toY) {
  const overlay = container.querySelector('#hover-overlay');
  const hLine   = container.querySelector('#h-line');
  const hDot    = container.querySelector('#h-dot');
  const tooltip = document.getElementById('svg-tooltip');
  if (!overlay || !tooltip) return;

  const svgEl = container.querySelector('svg');
  let ctm = null;

  overlay.addEventListener('mousemove', e => {
    if (!ctm) ctm = svgEl.getScreenCTM();
    const inv = ctm.inverse();
    const pt  = new DOMPoint(e.clientX, e.clientY).matrixTransform(inv);
    const rank = Math.max(0, Math.min(nG - 1,
      Math.round((pt.x - lm) / pw * (nG - 1))));

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
    tooltip.innerHTML =
      `<strong>${_esc(name)}</strong><br>Rank #${rank+1}<br>ES\u00a0${es.toFixed(4)}`;
  });

  overlay.addEventListener('mouseleave', () => {
    hLine.setAttribute('opacity', '0');
    hDot.setAttribute('opacity', '0');
    tooltip.style.display = 'none';
    ctm = null;
  });

  window.addEventListener('resize', () => { ctm = null; }, { passive: true });
}

// ── Plot header & stats ──────────────────────────────────────

export function updatePlotHeader(result) {
  document.getElementById('plot-name').textContent = result.name;
  document.getElementById('plot-section').style.display = 'block';
  const link = document.getElementById('plot-db-link');
  const url  = msigdbUrl(result.name, result.url);
  link.href  = url;
  link.style.display = 'inline-flex';
}

export function renderESStats(result, engine, showFDR) {
  const fmt = p => p == null ? '—'
    : Math.abs(p) < 0.001 ? p.toExponential(2) : p.toFixed(4);
  const sc = p => p != null && p < 0.05 ? 'sig' : '';

  const isPar = engine === 'parametric';
  const cells = [
    { l:'ES',     v:result.es.toFixed(4),  c:result.es>=0?'pos':'neg' },
    { l:'NES',    v:result.nes.toFixed(3), c:result.nes>=0?'pos':'neg' },
    { l:'NES-AD', v:result.nes_ad.toFixed(3), c:'' },
    { l:isPar?'p<sub>KS</sub>&thinsp;(Γ)':'p<sub>KS</sub>',
      v:fmt(result.pKS), c:sc(result.pKS) },
    { l:isPar?'p<sub>AD</sub>&thinsp;(GΓ)':'p<sub>AD</sub>',
      v:fmt(result.pAD), c:sc(result.pAD) },
    { l:'p<sub>Cauchy</sub>', v:fmt(result.pCauchy), c:sc(result.pCauchy) },
    ...(showFDR && result.fdr!=null
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
