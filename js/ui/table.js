// ═══════════════════════════════════════════════════════════
//  ui/table.js  ·  Results table
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

// Column definitions — order is render order
const COLS = [
  { id:'rank',    hdr:'#',                      extra:false, sort:false },
  { id:'name',    hdr:'Pathway',                extra:false, sort:true  },
  { id:'size',    hdr:'Size',                   extra:false, sort:true  },
  { id:'nes',     hdr:'NES',                    extra:false, sort:true  },
  { id:'nes_ad',  hdr:'NES-AD',                 extra:false, sort:true  },
  { id:'pKS',     hdr:'p<sub>KS</sub>',         extra:false, sort:true  },
  { id:'pAD',     hdr:'p<sub>AD</sub>',         extra:false, sort:true  },
  { id:'pCauchy', hdr:'p<sub>Cauchy</sub>',     extra:false, sort:true  },
  { id:'fdr',     hdr:'FDR',                    extra:false, sort:true, fdrOnly:true },
  // Extended
  { id:'es',      hdr:'ES',                     extra:true,  sort:true  },
  { id:'ad',      hdr:'AD',                     extra:true,  sort:true  },
  { id:'pKS_emp', hdr:'p<sub>KS</sub> (emp)',   extra:true,  sort:true  },
  { id:'pAD_emp', hdr:'p<sub>AD</sub> (emp)',   extra:true,  sort:true  },
];

// Sort state
let _sortCol = 'pCauchy';
let _sortDir = 1;           // 1 = ascending, -1 = descending

/**
 * Render (or re-render) the results table.
 *
 * @param {object[]} results   GSEA result array
 * @param {boolean}  showFDR  show FDR column (true when nPathways > 10)
 * @param {string}   engine   'gg' | 'empirical'
 * @param {Function} onSelect callback(result) when a row's plot button is clicked
 */
export function renderTable(results, showFDR, engine, onSelect) {
  _buildHead(showFDR, engine);
  _buildBody(results, showFDR, onSelect);
  _attachSort(results, showFDR, onSelect);
}

// ── Header ───────────────────────────────────────────────────
function _buildHead(showFDR, engine) {
  const isGG = engine === 'gg';
  const labelOverride = {
    pKS:    isGG ? 'p<sub>KS</sub>&thinsp;<small>(Γ)</small>'   : 'p<sub>KS</sub>',
    pAD:    isGG ? 'p<sub>AD</sub>&thinsp;<small>(GΓ)</small>'  : 'p<sub>AD</sub>',
  };

  const visibleCols = COLS.filter(c => {
    if (c.fdrOnly && !showFDR) return false;
    return true;
  });

  const ths = visibleCols.map(c => {
    const lbl   = labelOverride[c.id] ?? c.hdr;
    const xCls  = c.extra ? 'col-ext' : '';
    const sCls  = c.sort  ? 'sortable' : '';
    const arrow = c.sort  ? '<span class="si"></span>' : '';
    return `<th class="${xCls} ${sCls}" data-col="${c.id}">${lbl}${arrow}</th>`;
  });

  document.getElementById('rt-head').innerHTML = `<tr>${ths.join('')}</tr>`;
}

// ── Body ─────────────────────────────────────────────────────
function _buildBody(results, showFDR, onSelect) {
  const sorted = _sorted(results);
  const tbody  = document.getElementById('tbody');
  tbody.innerHTML = '';

  sorted.forEach((r, rank) => {
    const tr = document.createElement('tr');
    tr.dataset.name = r.name;

    const url   = msigdbUrl(r.name, r.url);
    const cells = _cells(r, rank, url, showFDR);
    tr.innerHTML = cells;

    // Plot button
    tr.querySelector('.btn-plot')?.addEventListener('click', e => {
      e.stopPropagation();
      _selectRow(tr);
      onSelect(r);
    });

    tbody.appendChild(tr);
  });

  // Auto-select + plot first row
  const first = tbody.firstElementChild;
  if (first) {
    first.classList.add('sel');
    const name = first.dataset.name;
    const r    = results.find(x => x.name === name);
    if (r) onSelect(r);
  }
}

function _cells(r, rank, url, showFDR) {
  const pc = p =>
    p == null ? '' :
    p < 0.001 ? 's001' :
    p < 0.01  ? 's01'  :
    p < 0.05  ? 's05'  : '';

  const fp = p =>
    p == null ? '—' :
    Math.abs(p) < 0.001 ? p.toExponential(2) : p.toFixed(4);

  // Build a map of cell HTML per column id
  const C = {
    rank:    `<td class="cn">${rank + 1}</td>`,
    name:    `<td class="cpname">
                <a class="path-link" href="${url}" target="_blank"
                   rel="noopener" title="${_esc(r.name)}">${_esc(r.name)}</a
                ><button class="btn-plot" title="View enrichment plot">📈</button>
              </td>`,
    size:    `<td class="num neu">${r.size}</td>`,
    nes:     `<td class="num ${r.nes    >= 0 ? 'pos' : 'neg'}">${r.nes.toFixed(3)}</td>`,
    nes_ad:  `<td class="num ${r.nes_ad >= 0 ? 'pos' : 'neg'}">${r.nes_ad.toFixed(3)}</td>`,
    pKS:     `<td class="pv ${pc(r.pKS)}">${fp(r.pKS)}</td>`,
    pAD:     `<td class="pv ${pc(r.pAD)}">${fp(r.pAD)}</td>`,
    pCauchy: `<td class="pv ${pc(r.pCauchy)}">${fp(r.pCauchy)}</td>`,
    fdr:     showFDR
               ? `<td class="fv ${r.fdr != null && r.fdr < 0.05 ? 'fsig' : ''}">${fp(r.fdr)}</td>`
               : '',
    // Extended
    es:      `<td class="num col-ext ${r.es >= 0 ? 'pos' : 'neg'}">${r.es.toFixed(4)}</td>`,
    ad:      `<td class="num col-ext neu">${r.ad.toFixed(2)}</td>`,
    pKS_emp: `<td class="pv col-ext ${pc(r.pKS_emp)}">${fp(r.pKS_emp)}</td>`,
    pAD_emp: `<td class="pv col-ext ${pc(r.pAD_emp)}">${fp(r.pAD_emp)}</td>`,
  };

  return COLS
    .filter(c => !(c.fdrOnly && !showFDR))
    .map(c => C[c.id] ?? '')
    .join('');
}

// ── Sort ─────────────────────────────────────────────────────
function _attachSort(results, showFDR, onSelect) {
  document.querySelectorAll('#rt-head th.sortable').forEach(th => {
    th.addEventListener('click', () => {
      const col = th.dataset.col;
      _sortDir = col === _sortCol ? -_sortDir : 1;
      _sortCol = col;
      // Update header arrows
      document.querySelectorAll('#rt-head th').forEach(h => {
        h.classList.remove('sort-asc', 'sort-desc');
      });
      th.classList.add(_sortDir === 1 ? 'sort-asc' : 'sort-desc');
      _buildBody(results, showFDR, onSelect);
    });
  });
}

function _sorted(results) {
  const key = {
    name:    r => r.name,
    size:    r => r.size,
    nes:     r => r.nes,
    nes_ad:  r => r.nes_ad,
    pKS:     r => r.pKS,
    pAD:     r => r.pAD,
    pCauchy: r => r.pCauchy,
    fdr:     r => r.fdr ?? 1,
    es:      r => r.es,
    ad:      r => r.ad,
    pKS_emp: r => r.pKS_emp,
    pAD_emp: r => r.pAD_emp,
  }[_sortCol] ?? (r => r.pCauchy);

  return [...results].sort((a, b) => {
    const av = key(a), bv = key(b);
    return typeof av === 'string'
      ? _sortDir * av.localeCompare(bv)
      : _sortDir * (av - bv);
  });
}

function _selectRow(tr) {
  document.querySelectorAll('#tbody tr').forEach(t => t.classList.remove('sel'));
  tr.classList.add('sel');
}

const _esc = s => (s ?? '').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
