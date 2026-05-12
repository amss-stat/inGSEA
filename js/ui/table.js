// ═══════════════════════════════════════════════════════════
//  ui/table.js  ·  v2.9 (Standard GSEA FDR Edition)
// ═══════════════════════════════════════════════════════════
'use strict';

import { msigdbUrl } from './fileio.js';

const COLS = [
  { id:'rank',     hdr:'#',                            extra:false, sort:false },
  { id:'name',     hdr:'Pathway',                      extra:false, sort:true  },
  { id:'size',     hdr:'Size',                         extra:false, sort:true  },
  // KS 核心指标
  { id:'nes',      hdr:'NES',                          extra:false, sort:true  },
  { id:'fdr_ks',   hdr:'FDR (KS)',                     extra:false, sort:true  },
  // AD 核心指标
  { id:'nes_ad',   hdr:'NES-AD',                       extra:false, sort:true  },
  { id:'fdr_ad',   hdr:'FDR (AD)',                     extra:false, sort:true  },
  // 综合指标
  { id:'pCauchy',  hdr:'p<sub>Cauchy</sub>',           extra:false, sort:true  },
  
  // 扩展列 (默认隐藏，点击 Extended columns 切换)
  { id:'pKS',      hdr:'p<sub>KS</sub>',               extra:true,  sort:true  },
  { id:'pAD',      hdr:'p<sub>AD</sub>',               extra:true,  sort:true  },
  { id:'pKS_emp',  hdr:'p<sub>KS</sub>&nbsp;(perm)',   extra:true,  sort:true  },
  { id:'pAD_emp',  hdr:'p<sub>AD</sub>&nbsp;(perm)',   extra:true,  sort:true  },
  { id:'pAD_par',  hdr:'p<sub>AD</sub>&nbsp;(par)',    extra:true,  sort:true  },
];

let _sortCol = 'fdr_ks'; // 默认按 KS 的 FDR 排序
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
    const cur = c.id === _sortCol ? (_sortDir === 1 ? 'sort-asc' : 'sort-desc') : '';
    return `<th class="${xc} ${sc} ${cur}" data-col="${c.id}">${lbl}${ar}</th>`;
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

    // 点击行选择通路
    tr.addEventListener('click', () => {
      _selectRow(tr); 
      onSelect(r);
    });

    tbody.appendChild(tr);
  });

  // 默认选中第一行
  const first = tbody.firstElementChild;
  if (first) {
    first.classList.add('sel');
    const r = results.find(x => x.name === first.dataset.name);
    if (r) onSelect(r);
  }
}

// ── Single cell renderer ──────────────────────────────────────
function _cell(id, r, rank, url) {
  // P-value 着色逻辑
  const pc = p =>
    p == null ? '' :
    p < 0.001 ? 's001' :
    p < 0.01  ? 's01'  :
    p < 0.05  ? 's05'  : '';

  // FDR 着色逻辑 (GSEA 标准: < 0.05 显著, < 0.25 值得关注)
  const fc = f =>
    f == null ? '' :
    f < 0.05  ? 'fsig' : 
    f < 0.25  ? 'fwarn' : '';

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

    case 'fdr_ks':
      return `<td class="fv ${fc(r.fdr_ks)}">${fp(r.fdr_ks)}</td>`;

    case 'nes_ad':
      return `<td class="num neu">${r.nes_ad.toFixed(3)}</td>`;

    case 'fdr_ad':
      return `<td class="fv ${fc(r.fdr_ad)}">${fp(r.fdr_ad)}</td>`;

    case 'pCauchy':
      return `<td class="pv ${pc(r.pCauchy)}">${fp(r.pCauchy)}</td>`;

    case 'pKS':
      return `<td class="pv col-ext ${pc(r.pKS)}">${fp(r.pKS)}</td>`;

    case 'pAD':
      return `<td class="pv col-ext ${pc(r.pAD)}">${fp(r.pAD)}</td>`;

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
      _buildHead(); // 更新箭头
      _buildBody(results, onSelect);
    });
  });
}

function _sorted(results) {
  const key = {
    name:    r => r.name,
    size:    r => r.size,
    nes:     r => r.nes,
    fdr_ks:  r => r.fdr_ks ?? 1,
    nes_ad:  r => r.nes_ad,
    fdr_ad:  r => r.fdr_ad ?? 1,
    pCauchy: r => r.pCauchy,
    pKS:     r => r.pKS,
    pAD:     r => r.pAD,
    pKS_emp: r => r.pKS_emp,
    pAD_emp: r => r.pAD_emp,
    pAD_par: r => r.pAD_par ?? 1,
  }[_sortCol] ?? (r => r.fdr_ks);

  return [...results].sort((a, b) => {
    const av = key(a), bv = key(b);
    if (typeof av === 'string') {
      return _sortDir * av.localeCompare(bv);
    }
    return _sortDir * (av - bv);
  });
}

function _selectRow(tr) {
  document.querySelectorAll('#tbody tr').forEach(t => t.classList.remove('sel'));
  tr.classList.add('sel');
}

const _esc = s =>
  (s ?? '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
