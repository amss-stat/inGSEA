// ═══════════════════════════════════════════════════════════
//  ui/table.js  ·  v0.3.0
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

const COLS = [
  { id:'rank',     hdr:'#',                            sort:false },
  { id:'name',     hdr:'Pathway',                      sort:true  },
  { id:'size',     hdr:'Size',                         sort:true  },
  { id:'nes',      hdr:'NES',                          sort:true  },
  { id:'nes_ad',   hdr:'NES-AD',                       sort:true  },
  { id:'pKS',      hdr:'p<sub>KS</sub>',               sort:true  },
  { id:'pAD',      hdr:'p<sub>AD</sub>',               sort:true  }, // Dynamic
  { id:'pCauchy',  hdr:'<strong>p<sub>Cauchy</sub></strong>', sort:true  },
  { id:'fdr_ks',   hdr:'FDR<sub>KS</sub>',             sort:true  },
  { id:'fdr_ad',   hdr:'<strong>FDR<sub>AD</sub></strong>', sort:true  }
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
    pAD: isPar
      ? '<strong>p<sub>AD</sub>&thinsp;<small>(par)</small></strong>'
      : '<strong>p<sub>AD</sub>&thinsp;<small>(perm)</small></strong>'
  };

  const ths = COLS.map(c => {
    const lbl = over[c.id] ?? c.hdr;
    const sc  = c.sort   ? 'sortable' : '';
    const ar  = c.sort   ? '<span class="si"></span>' : '';
    return `<th class="${sc}" data-col="${c.id}">${lbl}${ar}</th>`;
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

  const fc = f => 
    f == null ? '' :
    f < 0.05  ? 'fsig' : 
    f < 0.25  ? 'fwarn' : '';

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

    case 'nes_ad':
      return `<td class="num neu">${r.nes_ad.toFixed(3)}</td>`;

    case 'pKS':
      return `<td class="pv ${pc(r.pKS)}">${fp(r.pKS)}</td>`;

    case 'pAD':
      return `<td class="pv ${pc(r.pAD)}">${fp(r.pAD)}</td>`;

    case 'pCauchy':
      return `<td class="pv ${pc(r.pCauchy)}">${fp(r.pCauchy)}</td>`;

    case 'fdr_ks':
      return `<td class="fv ${fc(r.fdr_ks)}">${fp(r.fdr_ks)}</td>`;

    case 'fdr_ad':
      return `<td class="fv ${fc(r.fdr_ad)}">${fp(r.fdr_ad)}</td>`;

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
    fdr_ks:  r => r.fdr_ks ?? 1,
    fdr_ad:  r => r.fdr_ad ?? 1,
  }[_sortCol] ?? (r => r.pCauchy);

  return [...results].sort((a, b) => {
    const av = key(a), bv = key(b);
    
    // Safety net for null/undefined items
    if (av == null && bv != null) return 1;
    if (bv == null && av != null) return -1;
    if (av == null && bv == null) return 0;

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
