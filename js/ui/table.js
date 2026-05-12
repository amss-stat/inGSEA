// ═══════════════════════════════════════════════════════════
//  ui/table.js  ·  v2.6
//  Fixes:
//  • FDR col always rendered (col-ext); fp() never receives null
//    when results.length >= 2 (gsea.js now always sets fdr)
//  • showFDR flag controls FDR note in main.js only;
//    table always has the column available in extended view
//  • pKS_par / pAD_par shown in extended columns when available
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

const COLS = [
  { id:'rank',     hdr:'#',                            extra:false, sort:false },
  { id:'name',     hdr:'Pathway',                      extra:false, sort:true  },
  { id:'size',     hdr:'Size',                         extra:false, sort:true  },
  { id:'nes',      hdr:'NES',                          extra:false, sort:true  },
  { id:'nes_ad',   hdr:'NES-AD',                       extra:false, sort:true  },
  { id:'pKS',      hdr:'p<sub>KS</sub>',               extra:false, sort:true  },
  { id:'pAD',      hdr:'p<sub>AD</sub>',               extra:false, sort:true  },
  { id:'pCauchy',  hdr:'p<sub>Cauchy</sub>',           extra:false, sort:true  },
  // Extended columns
  { id:'fdr',      hdr:'FDR (BH)',                     extra:true,  sort:true  },
  { id:'es',       hdr:'ES',                           extra:true,  sort:true  },
  { id:'ad',       hdr:'AD',                           extra:true,  sort:true  },
  { id:'pKS_emp',  hdr:'p<sub>KS</sub>&nbsp;(perm)',   extra:true,  sort:true  },
  { id:'pAD_emp',  hdr:'p<sub>AD</sub>&nbsp;(perm)',   extra:true,  sort:true  },
  { id:'pKS_par',  hdr:'p<sub>KS</sub>&nbsp;(par)',    extra:true,  sort:true  },
  { id:'pAD_par',  hdr:'p<sub>AD</sub>&nbsp;(par)',    extra:true,  sort:true  },
];

let _sortCol = 'pCauchy';
let _sortDir = 1;

export function renderTable(results, showFDR, engine, onSelect) {
  _buildHead(engine);
  _buildBody(results, onSelect);
  _attachSort(results, onSelect);
}

// ── Head ──────────────────────────────────────────────────────
function _buildHead(engine) {
  const isPar = engine === 'parametric';
  const over  = {
    pKS: isPar
      ? 'p<sub>KS</sub>&thinsp;<small>(par)</small>'
      : 'p<sub>KS</sub>',
    pAD: isPar
      ? 'p<sub>AD</sub>&thinsp;<small>(par)</small>'
      : 'p<sub>AD</sub>',
  };

  const ths = COLS.map(c => {
    const lbl = over[c.id] ?? c.hdr;
    const xc  = c.extra  ? 'col-ext' : '';
    const sc  = c.sort   ? 'sortable' : '';
    const ar  = c.sort   ? '<span class="si"></span>' : '';
    return `<th class="${xc} ${sc}" data-col="${c.id}">${lbl}${ar}</th>`;
  });

  document.getElementById('rt-head').innerHTML =
    `<tr>${ths.join('')}</tr>`;
}

// ── Body ──────────────────────────────────────────────────────
function _buildBody(results, onSelect) {
  const sorted = _sorted(results);
  const tbody  = document.getElementById('tbody');
  tbody.innerHTML = '';

  sorted.forEach((r, rank) => {
    const tr = document.createElement('tr');
    tr.dataset.name = r.name;
    const url = msigdbUrl(r.name, r.url);
    tr.innerHTML = COLS.map(c => _cell(c.id, r, rank, url)).join('');

    tr.querySelector('.path-name-btn')?.addEventListener('click', e => {
      e.stopPropagation(); _selectRow(tr); onSelect(r);
    });
    tr.querySelector('.btn-plot')?.addEventListener('click', e => {
      e.stopPropagation(); _selectRow(tr); onSelect(r);
    });

    tbody.appendChild(tr);
  });

  const first = tbody.firstElementChild;
  if (first) {
    first.classList.add('sel');
    const r = results.find(x => x.name === first.dataset.name);
    if (r) onSelect(r);
  }
}

// ── Single cell renderer ──────────────────────────────────────
function _cell(id, r, rank, url) {
  const pc = p =>
    p == null ? '' :
    p < 0.001 ? 's001' :
    p < 0.01  ? 's01'  :
    p < 0.05  ? 's05'  : '';

  const fp = p =>
    p == null ? '—' :
    Math.abs(p) < 0.001 ? p.toExponential(2) : p.toFixed(4);

  switch (id) {
    case 'rank':
      return `<td class="cn">${rank + 1}</td>`;

    case 'name':
      return `<td class="cpname">
        <span class="path-name-btn" title="${_esc(r.name)}">${_esc(r.name)}</span>
        <button class="btn-plot" title="View enrichment curve">📈</button>
        <a class="btn-msigdb" href="${url}" target="_blank" rel="noopener"
           title="Jump to pathway database information">↗MSigDB</a>
      </td>`;

    case 'size':
      return `<td class="num neu">${r.size}</td>`;

    case 'nes':
      return `<td class="num ${r.nes >= 0 ? 'pos' : 'neg'}">${r.nes.toFixed(3)}</td>`;

    // Issue (5): NES-AD always neutral
    case 'nes_ad':
      return `<td class="num neu">${r.nes_ad.toFixed(3)}</td>`;

    case 'pKS':
      return `<td class="pv ${pc(r.pKS)}">${fp(r.pKS)}</td>`;

    case 'pAD':
      return `<td class="pv ${pc(r.pAD)}">${fp(r.pAD)}</td>`;

    case 'pCauchy':
      return `<td class="pv ${pc(r.pCauchy)}">${fp(r.pCauchy)}</td>`;

    case 'fdr':
      return `<td class="fv col-ext ${r.fdr != null && r.fdr < 0.05 ? 'fsig' : ''}">
        ${fp(r.fdr)}
      </td>`;

    case 'es':
      return `<td class="num col-ext ${r.es >= 0 ? 'pos' : 'neg'}">${r.es.toFixed(4)}</td>`;

    case 'ad':
      return `<td class="num col-ext neu">${r.ad.toFixed(2)}</td>`;

    case 'pKS_emp':
      return `<td class="pv col-ext ${pc(r.pKS_emp)}">${fp(r.pKS_emp)}</td>`;

    case 'pAD_emp':
      return `<td class="pv col-ext ${pc(r.pAD_emp)}">${fp(r.pAD_emp)}</td>`;

    case 'pKS_par':
      return `<td class="pv col-ext ${pc(r.pKS_par)}">${fp(r.pKS_par)}</td>`;

    case 'pAD_par':
      return `<td class="pv col-ext ${pc(r.pAD_par)}">${fp(r.pAD_par)}</td>`;

    default:
      return '<td></td>';
  }
}

// ── Sort ──────────────────────────────────────────────────────
function _attachSort(results, onSelect) {
  document.querySelectorAll('#rt-head th.sortable').forEach(th => {
    th.addEventListener('click', () => {
      const col = th.dataset.col;
      _sortDir  = col === _sortCol ? -_sortDir : 1;
      _sortCol  = col;
      document.querySelectorAll('#rt-head th')
        .forEach(h => h.classList.remove('sort-asc', 'sort-desc'));
      th.classList.add(_sortDir === 1 ? 'sort-asc' : 'sort-desc');
      _buildBody(results, onSelect);
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
    pKS_par: r => r.pKS_par ?? 1,
    pAD_par: r => r.pAD_par ?? 1,
  }[_sortCol] ?? (r => r.pCauchy);

  return [...results].sort((a, b) => {
    const av = key(a), bv = key(b);
    return typeof av === 'string'
      ? _sortDir * av.localeCompare(bv)
      : _sortDir * (av - bv);
  });
}

function _selectRow(tr) {
  document.querySelectorAll('#tbody tr')
    .forEach(t => t.classList.remove('sel'));
  tr.classList.add('sel');
}

const _esc = s =>
  (s ?? '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
