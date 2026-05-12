// ═══════════════════════════════════════════════════════════
//  ui/table.js  ·  Results table renderer & sorter
// ═══════════════════════════════════════════════════════════
'use strict';

import { inferMSigDBUrl } from './fileio.js';

const DEFAULT_COLS = [
  { id: 'rank',     label: '#',              extra: false, sortable: false },
  { id: 'name',     label: 'Pathway',        extra: false, sortable: true  },
  { id: 'size',     label: 'Size',           extra: false, sortable: true  },
  { id: 'nes',      label: 'NES',            extra: false, sortable: true  },
  { id: 'nes_ad',   label: 'NES-AD',         extra: false, sortable: true  },
  { id: 'pKS',      label: 'p<sub>KS</sub>', extra: false, sortable: true  },
  { id: 'pAD',      label: 'p<sub>AD</sub>', extra: false, sortable: true  },
  { id: 'pCauchy',  label: 'p<sub>Cauchy</sub>', extra: false, sortable: true },
  { id: 'fdr',      label: 'FDR',            extra: false, sortable: true, fdrCol: true },
  // Extended (hidden by default)
  { id: 'es',       label: 'ES',             extra: true,  sortable: true  },
  { id: 'ad',       label: 'AD',             extra: true,  sortable: true  },
  { id: 'pKS_emp',  label: 'p<sub>KS</sub> (emp)', extra: true, sortable: true },
  { id: 'pAD_emp',  label: 'p<sub>AD</sub> (emp)', extra: true, sortable: true },
];

let sortCol = 'pCauchy', sortDir = 1;   // 1=asc, -1=desc

/**
 * Render the full results table.
 *
 * @param {object[]} results
 * @param {boolean}  showFDR    – if >10 pathways
 * @param {string}   engine     – 'gg' | 'empirical'
 * @param {Function} onSelect   – callback(result) when pathway name clicked
 */
export function renderTable(results, showFDR, engine, onSelect) {
  _buildHeader(showFDR, engine);
  _buildBody(results, showFDR, engine, onSelect);
  _attachSort(results, showFDR, engine, onSelect);
}

function _buildHeader(showFDR, engine) {
  const thead = document.getElementById('rt-head');
  const visibleCols = DEFAULT_COLS.filter(c => {
    if (c.fdrCol && !showFDR) return false;
    return true;
  });

  // Column label overrides for engine
  const labelMap = {
    pKS:     engine === 'gg' ? 'p<sub>KS</sub> <span style="font-size:8px;color:var(--muted)">(Γ)</span>'     : 'p<sub>KS</sub>',
    pAD:     engine === 'gg' ? 'p<sub>AD</sub> <span style="font-size:8px;color:var(--muted)">(GΓ)</span>'    : 'p<sub>AD</sub>',
    pCauchy: engine === 'gg' ? 'p<sub>Cauchy</sub>' : 'p<sub>Cauchy</sub>',
  };

  thead.innerHTML = `<tr>${visibleCols.map(c => {
    const lbl   = labelMap[c.id] || c.label;
    const extra = c.extra ? 'col-extra' : '';
    const sort  = c.sortable ? 'sortable' : '';
    const arrow = c.sortable ? '<span class="sort-icon"></span>' : '';
    return `<th class="${extra} ${sort}" data-col="${c.id}">${lbl}${arrow}</th>`;
  }).join('')}</tr>`;
}

function _buildBody(results, showFDR, engine, onSelect) {
  const tbody  = document.getElementById('tbody');
  const sorted = _sortResults([...results], sortCol, sortDir);
  tbody.innerHTML = '';

  sorted.forEach((r, i) => {
    const tr = document.createElement('tr');
    tr.dataset.name = r.name;

    const url     = inferMSigDBUrl(r.name, r.url);
    const pKS     = r.pKS_fit,  pAD = r.pAD_fit,  pC = r.pCauchy;
    const pKSe    = r.pKS_emp,  pADe = r.pAD_emp;
    const fdr     = r.fdr;

    const pcls = p => p < 0.001 ? 'sig001' : p < 0.01 ? 'sig01' : p < 0.05 ? 'sig05' : '';

    // Build cells in column order
    const cells = {
      rank:    `<td style="color:var(--muted);font-size:11px">${i + 1}</td>`,
      name:    `<td class="pname">
                  <a class="path-link" href="${url}" target="_blank"
                     title="${esc(r.name)}">${esc(r.name)}</a>
                  <span class="view-curve" data-name="${esc(r.name)}" title="View enrichment plot">📈</span>
                </td>`,
      size:    `<td class="num" style="color:var(--text2)">${r.size}</td>`,
      nes:     `<td class="num ${r.nes  >= 0 ? 'pos' : 'neg'}">${r.nes.toFixed(3)}</td>`,
      nes_ad:  `<td class="num ${r.nes_ad >= 0 ? 'pos' : 'neg'}">${r.nes_ad.toFixed(3)}</td>`,
      pKS:     `<td class="pval ${pcls(pKS)}">${fmtP(pKS)}</td>`,
      pAD:     `<td class="pval ${pcls(pAD)}">${fmtP(pAD)}</td>`,
      pCauchy: `<td class="pval ${pcls(pC)}">${fmtP(pC)}</td>`,
      fdr:     showFDR ? `<td class="fdr-val ${fdr != null && fdr < 0.05 ? 'sig' : ''}">${fdr != null ? fmtP(fdr) : '—'}</td>` : '',
      // Extended
      es:      `<td class="num col-extra ${r.es >= 0 ? 'pos' : 'neg'}">${r.es.toFixed(4)}</td>`,
      ad:      `<td class="num col-extra">${r.ad.toFixed(2)}</td>`,
      pKS_emp: `<td class="pval col-extra ${pcls(pKSe)}">${fmtP(pKSe)}</td>`,
      pAD_emp: `<td class="pval col-extra ${pcls(pADe)}">${fmtP(pADe)}</td>`,
    };

    const visibleCols = DEFAULT_COLS.filter(c => !(c.fdrCol && !showFDR));
    tr.innerHTML = visibleCols.map(c => cells[c.id] || '').join('');

    // "View curve" click
    tr.querySelector('.view-curve')?.addEventListener('click', e => {
      e.preventDefault();
      e.stopPropagation();
      _selectRow(tr);
      onSelect(r);
    });

    tbody.appendChild(tr);
  });

  // Auto-select first row
  const first = tbody.firstElementChild;
  if (first) {
    first.classList.add('sel');
    const firstName = first.dataset.name;
    const firstResult = results.find(r => r.name === firstName);
    if (firstResult) onSelect(firstResult);
  }
}

function _selectRow(tr) {
  document.querySelectorAll('#tbody tr').forEach(t => t.classList.remove('sel'));
  tr.classList.add('sel');
}

function _attachSort(results, showFDR, engine, onSelect) {
  document.querySelectorAll('#rt-head th.sortable').forEach(th => {
    th.addEventListener('click', () => {
      const col = th.dataset.col;
      if (sortCol === col) {
        sortDir = -sortDir;
      } else {
        sortCol = col;
        sortDir = 1;
      }
      // Update header arrows
      document.querySelectorAll('#rt-head th').forEach(h => {
        h.classList.remove('sort-asc', 'sort-desc');
      });
      th.classList.add(sortDir === 1 ? 'sort-asc' : 'sort-desc');
      _buildBody(results, showFDR, engine, onSelect);
    });
  });
}

function _sortResults(res, col, dir) {
  const key = {
    name:    r => r.name,
    size:    r => r.size,
    nes:     r => r.nes,
    nes_ad:  r => r.nes_ad,
    pKS:     r => r.pKS_fit,
    pAD:     r => r.pAD_fit,
    pCauchy: r => r.pCauchy,
    fdr:     r => r.fdr ?? 1,
    es:      r => r.es,
    ad:      r => r.ad,
    pKS_emp: r => r.pKS_emp,
    pAD_emp: r => r.pAD_emp,
  }[col] || (r => r.pCauchy);

  return res.sort((a, b) => {
    const av = key(a), bv = key(b);
    if (typeof av === 'string') return dir * av.localeCompare(bv);
    return dir * (av - bv);
  });
}

export function fmtP(p) {
  if (p == null) return '—';
  return Math.abs(p) < 0.001 ? p.toExponential(2) : p.toFixed(4);
}

function esc(s) {
  return (s || '').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}
