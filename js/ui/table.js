// ═══════════════════════════════════════════════════════════
//  ui/table.js  ·  v2.9 (GSEA standard FDR Edition)
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

// 定义表格列：移除 ES 和 AD，增加 fdr_ks 和 fdr_ad
const COLS = [
  { id:'rank',     hdr:'#',                            extra:false, sort:false },
  { id:'name',     hdr:'Pathway',                      extra:false, sort:true  },
  { id:'size',     hdr:'Size',                         extra:false, sort:true  },
  { id:'nes',      hdr:'NES',                          extra:false, sort:true  },
  { id:'fdr_ks',   hdr:'FDR (KS)',                     extra:false, sort:true  },
  { id:'nes_ad',   hdr:'NES-AD',                       extra:false, sort:true  },
  { id:'fdr_ad',   hdr:'FDR (AD)',                     extra:false, sort:true  },
  // 辅助/次要信息放入扩展列
  { id:'pKS',      hdr:'p<sub>KS</sub>',               extra:true,  sort:true  },
  { id:'pAD',      hdr:'p<sub>AD</sub>',               extra:true,  sort:true  },
  { id:'pCauchy',  hdr:'p<sub>Cauchy</sub>',           extra:true,  sort:true  },
  { id:'pKS_emp',  hdr:'p<sub>KS</sub>&nbsp;(perm)',   extra:true,  sort:true  },
  { id:'pAD_emp',  hdr:'p<sub>AD</sub>&nbsp;(perm)',   extra:true,  sort:true  },
  { id:'pAD_par',  hdr:'p<sub>AD</sub>&nbsp;(par)',    extra:true,  sort:true  },
];

let _sortCol = 'fdr_ks'; // 默认按标准 FDR 排序
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
    pKS: 'p<sub>KS</sub>&thinsp;<small>(perm)</small>',
    pAD: isPar ? 'p<sub>AD</sub>&thinsp;<small>(par)</small>' : 'p<sub>AD</sub>',
  };

  const ths = COLS.map(c => {
    const lbl = over[c.id] ?? c.hdr;
    const xc  = c.extra  ? 'col-ext' : '';
    const sc  = c.sort   ? 'sortable' : '';
    const ar  = c.sort   ? '<span class="si"></span>' : '';
    return `<th class="${xc} ${sc}" data-col="${c.id}">${lbl}${ar}</th>`;
  });

  document.getElementById('rt-head').innerHTML = `<tr>${ths.join('')}</tr>`;
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
  // P-value 样式 (0.05/0.01)
  const pc = p =>
    p == null ? '' :
    p < 0.01  ? 's01'  :
    p < 0.05  ? 's05'  : '';

  // FDR 样式：GSEA 官方标准中 0.25 也是一个重要阈值
  const fc = f =>
    f == null ? '' :
    f < 0.05  ? 'fsig' :  // 强显著 (绿色/深色)
    f < 0.25  ? 'fwarn' : ''; // 潜在显著 (橙色/浅色)

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
           title="Jump to MSigDB">↗</a>
      </td>`;

    case 'size':
      return `<td class="num neu">${r.size}</td>`;

    case 'nes':
      return `<td class="num ${r.nes >= 0 ? 'pos' : 'neg'}">${r.nes.toFixed(3)}</td>`;

    case 'nes_ad':
      return `<td class="num neu">${r.nes_ad.toFixed(3)}</td>`;

    case 'fdr_ks':
      return `<td class="fv ${fc(r.fdr_ks)}">${fp(r.fdr_ks)}</td>`;

    case 'fdr_ad':
      return `<td class="fv ${fc(r.fdr_ad)}">${fp(r.fdr_ad)}</td>`;

    case 'pKS':
      return `<td class="pv ${pc(r.pKS)}">${fp(r.pKS)}</td>`;

    case 'pAD':
      return `<td class="pv ${pc(r.pAD)}">${fp(r.pAD)}</td>`;

    case 'pCauchy':
      return `<td class="pv col-ext ${pc(r.pCauchy)}">${fp(r.pCauchy)}</td>`;

    case 'pKS_emp':
      return `<td class="pv col-ext ${pc(r.pKS_emp)}">${fp(r.pKS_emp)}</td>`;

    case 'pAD_emp':
      return `<td class="pv col-ext ${pc(r.pAD_emp)}">${fp(r.pAD_emp)}</td>`;

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
    fdr_ks:  r => r.fdr_ks,
    fdr_ad:  r => r.fdr_ad,
    pKS:     r => r.pKS,
    pAD:     r => r.pAD,
    pCauchy: r => r.pCauchy,
    pKS_emp: r => r.pKS_emp,
    pAD_emp: r => r.pAD_emp,
    pAD_par: r => r.pAD_par ?? 1,
  }[_sortCol] ?? (r => r.fdr_ks);

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
