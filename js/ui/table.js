// ═══════════════════════════════════════════════════════════
//  ui/table.js  ·  Results table
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

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
  { id:'es',      hdr:'ES',                     extra:true,  sort:true  },
  { id:'ad',      hdr:'AD',                     extra:true,  sort:true  },
  { id:'pKS_emp', hdr:'p<sub>KS</sub> (perm)',  extra:true,  sort:true  },
  { id:'pAD_emp', hdr:'p<sub>AD</sub> (perm)',  extra:true,  sort:true  },
];

let _sortCol = 'pCauchy';
let _sortDir = 1;

export function renderTable(results, showFDR, engine, onSelect) {
  _buildHead(showFDR, engine);
  _buildBody(results, showFDR, onSelect);
  _attachSort(results, showFDR, onSelect);
}

function _buildHead(showFDR, engine) {
  const isPar = engine === 'parametric';
  const over = {
    pKS: isPar ? 'p<sub>KS</sub>&thinsp;<small>(Γ)</small>'  : 'p<sub>KS</sub>',
    pAD: isPar ? 'p<sub>AD</sub>&thinsp;<small>(GΓ)</small>' : 'p<sub>AD</sub>',
  };

  const vis = COLS.filter(c => !(c.fdrOnly && !showFDR));
  const ths = vis.map(c => {
    const lbl = over[c.id] ?? c.hdr;
    const xc  = c.extra ? 'col-ext' : '';
    const sc  = c.sort  ? 'sortable' : '';
    const ar  = c.sort  ? '<span class="si"></span>' : '';
    return `<th class="${xc} ${sc}" data-col="${c.id}">${lbl}${ar}</th>`;
  });

  document.getElementById('rt-head').innerHTML = `<tr>${ths.join('')}</tr>`;
}

function _buildBody(results, showFDR, onSelect) {
  const sorted = _sorted(results);
  const tbody  = document.getElementById('tbody');
  tbody.innerHTML = '';

  sorted.forEach((r, rank) => {
    const tr  = document.createElement('tr');
    tr.dataset.name = r.name;

    const url = msigdbUrl(r.name, r.url);
    const C   = _cells(r, rank, url, showFDR);
    const vis = COLS.filter(c => !(c.fdrOnly && !showFDR));
    tr.innerHTML = vis.map(c => C[c.id] ?? '').join('');

    tr.querySelector('.btn-plot')?.addEventListener('click', e => {
      e.stopPropagation();
      _selectRow(tr);
      onSelect(r);
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

function _cells(r, rank, url, showFDR) {
  const pc = p =>
    p == null ? '' :
    p < 0.001 ? 's001' :
    p < 0.01  ? 's01'  :
    p < 0.05  ? 's05'  : '';

  const fp = p =>
    p == null ? '—' :
    Math.abs(p) < 0.001 ? p.toExponential(2) : p.toFixed(4);

  return {
    rank:    `<td class="cn">${rank + 1}</td>`,
    name:    `<td class="cpname"><a class="path-link" href="${url}" target="_blank" rel="noopener" title="${_esc(r.name)}">${_esc(r.name)}</a><button class="btn-plot" title="View enrichment plot">📈</button></td>`,
    size:    `<td class="num neu">${r.size}</td>`,
    nes:     `<td class="num ${r.nes>=0?'pos':'neg'}">${r.nes.toFixed(3)}</td>`,
    nes_ad:  `<td class="num ${r.nes_ad>=0?'pos':'neg'}">${r.nes_ad.toFixed(3)}</td>`,
    pKS:     `<td class="pv ${pc(r.pKS)}">${fp(r.pKS)}</td>`,
    pAD:     `<td class="pv ${pc(r.pAD)}">${fp(r.pAD)}</td>`,
    pCauchy: `<td class="pv ${pc(r.pCauchy)}">${fp(r.pCauchy)}</td>`,
    fdr:     showFDR
               ? `<td class="fv ${r.fdr!=null&&r.fdr<0.05?'fsig':''}">${fp(r.fdr)}</td>`
               : '',
    es:      `<td class="num col-ext ${r.es>=0?'pos':'neg'}">${r.es.toFixed(4)}</td>`,
    ad:      `<td class="num col-ext neu">${r.ad.toFixed(2)}</td>`,
    pKS_emp: `<td class="pv col-ext ${pc(r.pKS_emp)}">${fp(r.pKS_emp)}</td>`,
    pAD_emp: `<td class="pv col-ext ${pc(r.pAD_emp)}">${fp(r.pAD_emp)}</td>`,
  };
}

function _attachSort(results, showFDR, onSelect) {
  document.querySelectorAll('#rt-head th.sortable').forEach(th => {
    th.addEventListener('click', () => {
      const col = th.dataset.col;
      _sortDir = col === _sortCol ? -_sortDir : 1;
      _sortCol = col;
      document.querySelectorAll('#rt-head th').forEach(
        h => h.classList.remove('sort-asc','sort-desc'));
      th.classList.add(_sortDir === 1 ? 'sort-asc' : 'sort-desc');
      _buildBody(results, showFDR, onSelect);
    });
  });
}

function _sorted(results) {
  const key = {
    name:r=>r.name, size:r=>r.size, nes:r=>r.nes, nes_ad:r=>r.nes_ad,
    pKS:r=>r.pKS, pAD:r=>r.pAD, pCauchy:r=>r.pCauchy,
    fdr:r=>r.fdr??1, es:r=>r.es, ad:r=>r.ad,
    pKS_emp:r=>r.pKS_emp, pAD_emp:r=>r.pAD_emp
  }[_sortCol] || (r=>r.pCauchy);
  return [...results].sort((a,b) => {
    const av=key(a), bv=key(b);
    return typeof av==='string' ? _sortDir*av.localeCompare(bv) : _sortDir*(av-bv);
  });
}

function _selectRow(tr) {
  document.querySelectorAll('#tbody tr').forEach(t=>t.classList.remove('sel'));
  tr.classList.add('sel');
}

const _esc = s => (s??'').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
