// ═══════════════════════════════════════════════════════════
//  ui/plot.js  ·  SVG enrichment-curve renderer
// ═══════════════════════════════════════════════════════════
'use strict';

import { inferMSigDBUrl } from './fileio.js';

const MARGIN = { left: 62, right: 18, top: 22, bottom: 50, hitH: 26 };
const SVG_H  = 320;

/**
 * Draw enrichment curve as an inline SVG.
 *
 * @param {object} result   – single GSEA result object
 * @param {object[]} pathways – full pathway list (for mask lookup)
 * @param {string[]} geneNames
 * @param {HTMLElement} container
 */
export function drawEnrichmentSVG(result, pathways, geneNames, container) {
  const W  = container.clientWidth || 720;
  const H  = SVG_H;
  const ml = MARGIN.left, mr = MARGIN.right;
  const mt = MARGIN.top,  mb = MARGIN.bottom;
  const hitH = MARGIN.hitH;
  const pw = W - ml - mr;
  const ph = H - mt - mb - hitH - 6;

  const curve   = result.curve;
  const nG      = curve.length;
  const ord     = result.obsOrd;
  const isPos   = result.es >= 0;
  const accentC = isPos ? '#1e7e4a' : '#b91c1c';
  const fillC   = isPos ? 'rgba(30,126,74,.12)' : 'rgba(185,28,28,.10)';

  // Y range
  let lo = 0, hi = 0;
  for (let i = 0; i < nG; i++) { if (curve[i] < lo) lo = curve[i]; if (curve[i] > hi) hi = curve[i]; }
  const pad  = Math.max(Math.abs(hi), Math.abs(lo)) * 0.15 + 0.02;
  const yMin = lo - pad, yMax = hi + pad, yR = yMax - yMin;

  const toX = i  => ml + (i / (nG - 1)) * pw;
  const toY = es => mt + ph * (1 - (es - yMin) / yR);
  const y0  = toY(0);

  // Pathway mask lookup
  const pathway = pathways.find(p => p.name === result.name);
  const mask    = pathway ? pathway.mask : null;

  // ── Build path data ──
  let pathD = `M ${toX(0).toFixed(2)},${toY(curve[0]).toFixed(2)}`;
  // Subsample for very large gene sets (>5000 genes) to keep SVG compact
  const step = nG > 5000 ? Math.ceil(nG / 5000) : 1;
  for (let i = step; i < nG; i += step) {
    pathD += ` L ${toX(i).toFixed(2)},${toY(curve[i]).toFixed(2)}`;
  }
  if (step > 1) pathD += ` L ${toX(nG - 1).toFixed(2)},${toY(curve[nG - 1]).toFixed(2)}`;

  // ── Fill path ──
  const fillD = `M ${toX(0).toFixed(2)},${y0.toFixed(2)} ` + pathD.slice(1)
    + ` L ${toX(nG - 1).toFixed(2)},${y0.toFixed(2)} Z`;

  // ── Peak marker ──
  const peakX = toX(result.peakIdx).toFixed(2);
  const peakY = toY(curve[result.peakIdx]).toFixed(2);

  // ── Y-axis ticks ──
  const nTicks = 5;
  const tickVals = Array.from({ length: nTicks + 1 }, (_, i) => yMin + yR * i / nTicks);

  // ── Gene-hit strips ── (only draw hits, no loop over non-hits)
  let hitStrips = '';
  if (mask) {
    const hy = mt + ph + 6;
    for (let i = 0; i < nG; i++) {
      if (mask[ord[i]]) {
        const x = toX(i).toFixed(2);
        hitStrips += `<line class="hit-strip" x1="${x}" y1="${hy}" x2="${x}" y2="${hy + hitH - 4}"
          stroke="${accentC}" stroke-width="1.2" opacity="0.7"/>`;
      }
    }
  }

  // ── Invisible hover overlay lines ── (for tooltip)
  // We use a thin transparent rect overlay and JS mouse-move for performance
  // (adding nGenes <rect> elements is too slow for large sets).

  // ── SVG markup ──
  const svg = `
<svg xmlns="http://www.w3.org/2000/svg"
     viewBox="0 0 ${W} ${H}"
     width="${W}" height="${H}"
     style="font-family:'Inter',sans-serif;user-select:none"
     id="es-svg">

  <defs>
    <clipPath id="plot-clip">
      <rect x="${ml}" y="${mt}" width="${pw}" height="${ph}"/>
    </clipPath>
  </defs>

  <!-- Background -->
  <rect x="0" y="0" width="${W}" height="${H}" fill="#fafbfc"/>
  <rect x="${ml}" y="${mt}" width="${pw}" height="${ph}" fill="white" stroke="#dde1e7" stroke-width="0.7"/>

  <!-- Grid lines -->
  ${tickVals.map(v => {
    const y = toY(v).toFixed(2);
    return `<line x1="${ml}" y1="${y}" x2="${ml + pw}" y2="${y}"
      stroke="#e8ecf0" stroke-width="0.6" stroke-dasharray="3,4"/>`;
  }).join('\n  ')}
  ${[1, 2, 3, 4].map(k => {
    const x = (ml + pw * k / 5).toFixed(2);
    return `<line x1="${x}" y1="${mt}" x2="${x}" y2="${mt + ph}"
      stroke="#e8ecf0" stroke-width="0.6" stroke-dasharray="3,4"/>`;
  }).join('\n  ')}

  <!-- Zero line -->
  <line x1="${ml}" y1="${y0.toFixed(2)}" x2="${ml + pw}" y2="${y0.toFixed(2)}"
        stroke="#c8cdd6" stroke-width="1"/>

  <!-- Fill -->
  <path d="${fillD}" fill="${fillC}" clip-path="url(#plot-clip)"/>

  <!-- ES curve -->
  <path d="${pathD}" fill="none" stroke="${accentC}" stroke-width="1.8"
        stroke-linejoin="round" stroke-linecap="round"
        clip-path="url(#plot-clip)" id="es-line"/>

  <!-- Peak marker -->
  <line x1="${peakX}" y1="${mt}" x2="${peakX}" y2="${mt + ph}"
        stroke="#b45309" stroke-width="1" stroke-dasharray="4,3" opacity="0.7"/>
  <circle cx="${peakX}" cy="${peakY}" r="4" fill="#b45309" stroke="white" stroke-width="1.2"/>

  <!-- Gene-hit strip area -->
  <rect x="${ml}" y="${mt + ph + 6}" width="${pw}" height="${hitH}"
        fill="#f0f2f5" stroke="#dde1e7" stroke-width="0.5"/>
  ${hitStrips}
  <text x="${ml + pw / 2}" y="${mt + ph + 6 + hitH - 5}"
        text-anchor="middle" fill="#9aa3b0" font-size="8">GENE HITS</text>

  <!-- Axes -->
  <line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}"
        stroke="#adb5bd" stroke-width="1.2"/>
  <line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}"
        stroke="#adb5bd" stroke-width="1.2"/>

  <!-- Y-axis labels -->
  ${tickVals.map(v => `
    <text x="${ml - 6}" y="${(toY(v) + 3.5).toFixed(2)}"
          text-anchor="end" fill="#6b7585" font-size="9.5">${v.toFixed(2)}</text>`
  ).join('')}

  <!-- Y-axis title -->
  <text transform="translate(13,${mt + ph / 2}) rotate(-90)"
        text-anchor="middle" fill="#6b7585" font-size="10">Enrichment Score</text>

  <!-- X-axis labels -->
  <text x="${ml}" y="${mt + ph + hitH + 20}" text-anchor="middle"
        fill="#6b7585" font-size="9.5">1</text>
  <text x="${(ml + pw / 2).toFixed(2)}" y="${mt + ph + hitH + 20}"
        text-anchor="middle" fill="#6b7585" font-size="9.5">${Math.floor(nG / 2)}</text>
  <text x="${ml + pw}" y="${mt + ph + hitH + 20}"
        text-anchor="middle" fill="#6b7585" font-size="9.5">${nG}</text>
  <text x="${(ml + pw / 2).toFixed(2)}" y="${mt + ph + hitH + 36}"
        text-anchor="middle" fill="#6b7585" font-size="11" font-weight="500">Gene rank</text>

  <!-- Hover overlay (transparent, for JS tooltip) -->
  <rect id="hover-overlay"
        x="${ml}" y="${mt}" width="${pw}" height="${ph}"
        fill="transparent" cursor="crosshair"/>

  <!-- Hover indicator line (hidden by default) -->
  <line id="hover-line" x1="0" y1="${mt}" x2="0" y2="${mt + ph}"
        stroke="#1a5fa8" stroke-width="1" stroke-dasharray="3,3"
        opacity="0" pointer-events="none"/>
  <circle id="hover-dot" cx="0" cy="0" r="3.5"
          fill="#1a5fa8" stroke="white" stroke-width="1"
          opacity="0" pointer-events="none"/>
</svg>`;

  container.innerHTML = svg;

  // ── Attach hover tooltip ──
  attachSVGTooltip(container, curve, ord, geneNames, nG, ml, pw, mt, ph, toX, toY);
}

/** Attach mouse-move tooltip to SVG hover overlay. */
function attachSVGTooltip(container, curve, ord, geneNames, nG, ml, pw, mt, ph, toX, toY) {
  const overlay = container.querySelector('#hover-overlay');
  const hLine   = container.querySelector('#hover-line');
  const hDot    = container.querySelector('#hover-dot');
  if (!overlay) return;

  let tooltip = document.getElementById('svg-tooltip');
  if (!tooltip) {
    tooltip = document.createElement('div');
    tooltip.id = 'svg-tooltip';
    tooltip.className = 'svg-tooltip';
    document.body.appendChild(tooltip);
  }

  // Get SVG bounding rect for coordinate conversion
  const getSVGRect = () => container.querySelector('svg').getBoundingClientRect();

  overlay.addEventListener('mousemove', e => {
    const rect  = getSVGRect();
    const svgEl = container.querySelector('svg');
    const vb    = svgEl.viewBox.baseVal;
    const scaleX = vb.width  / rect.width;
    const scaleY = vb.height / rect.height;

    const svgX = (e.clientX - rect.left) * scaleX;
    // const svgY = (e.clientY - rect.top)  * scaleY;

    // Map X to gene rank
    const rank = Math.max(0, Math.min(nG - 1, Math.round((svgX - ml) / pw * (nG - 1))));
    const geneIdx = ord[rank];
    const geneName = geneNames[geneIdx];
    const es = curve[rank];

    const lx = toX(rank).toFixed(2);
    const ly = toY(es).toFixed(2);

    hLine.setAttribute('x1', lx);
    hLine.setAttribute('x2', lx);
    hLine.setAttribute('opacity', '0.6');
    hDot.setAttribute('cx', lx);
    hDot.setAttribute('cy', ly);
    hDot.setAttribute('opacity', '1');

    tooltip.style.display = 'block';
    tooltip.style.left    = (e.clientX + 12) + 'px';
    tooltip.style.top     = (e.clientY - 36) + 'px';
    tooltip.innerHTML =
      `<strong>${geneName}</strong><br>Rank: #${rank + 1}<br>ES: ${es.toFixed(4)}`;
  });

  overlay.addEventListener('mouseleave', () => {
    hLine.setAttribute('opacity', '0');
    hDot.setAttribute('opacity', '0');
    tooltip.style.display = 'none';
  });
}

/** Update plot-panel header with pathway name & MSigDB link. */
export function updatePlotHeader(result) {
  document.getElementById('plot-name').textContent = result.name;

  const link = document.getElementById('plot-db-link');
  const url  = inferMSigDBUrl(result.name, result.url);
  if (url) {
    link.href = url;
    link.style.display = 'inline-flex';
  } else {
    link.style.display = 'none';
  }

  document.getElementById('plot-section').style.display = 'block';
}

/** Render the ES-stats row below the SVG. */
export function renderESStats(result, engine) {
  const fmt  = v => Math.abs(v) < 0.001 ? v.toExponential(2) : v.toFixed(4);
  const sigCls = p => p < 0.01 ? 'sig' : '';

  const pKS  = result.pKS_fit;
  const pAD  = result.pAD_fit;
  const pC   = result.pCauchy;
  const fdr  = result.fdr != null ? result.fdr : null;

  const cells = [
    { lbl: 'ES',      val: result.es.toFixed(4),  cls: result.es >= 0 ? 'pos' : 'neg' },
    { lbl: 'NES',     val: result.nes.toFixed(3),  cls: result.nes >= 0 ? 'pos' : 'neg' },
    { lbl: 'NES-AD',  val: result.nes_ad.toFixed(3), cls: '' },
    { lbl: engine === 'gg' ? 'p<sub>KS</sub> (Γ)' : 'p<sub>KS</sub>',
      val: fmt(pKS), cls: sigCls(pKS) },
    { lbl: engine === 'gg' ? 'p<sub>AD</sub> (GΓ)' : 'p<sub>AD</sub>',
      val: fmt(pAD), cls: sigCls(pAD) },
    { lbl: 'p<sub>Cauchy</sub>', val: fmt(pC), cls: sigCls(pC) },
    ...(fdr != null ? [{ lbl: 'FDR', val: fmt(fdr), cls: sigCls(fdr) }] : []),
    { lbl: 'Size',    val: result.size + ' genes', cls: '' }
  ];

  document.getElementById('es-stats-row').innerHTML = cells.map(c => `
    <div class="es-cell">
      <span class="es-lbl">${c.lbl}</span>
      <span class="es-val ${c.cls}">${c.val}</span>
    </div>`
  ).join('');
}
