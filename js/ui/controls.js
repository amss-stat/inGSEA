// ═══════════════════════════════════════════════════════════
//  ui/controls.js  ·  v0.3.0
// ═══════════════════════════════════════════════════════════
'use strict';

import { buildMasks } from './fileio.js';

// ── Log ──────────────────────────────────────────────────────
export function log(msg, cls = '') {
  const el = document.getElementById('log');
  const p  = document.createElement('p');
  if (cls) p.className = cls;
  const t = new Date().toLocaleTimeString('en', { hour12: false });
  p.textContent = `[${t}] ${msg}`;
  el.appendChild(p);
  el.scrollTop = el.scrollHeight;
}

// ── Progress ─────────────────────────────────────────────────
export function setProgress(pct, phase, msg) {
  document.getElementById('prog-fill').style.width =
    Math.min(pct, 100) + '%';
  document.getElementById('prog-pct').textContent  = pct + '%';
  if (phase) document.getElementById('prog-phase').textContent = phase;
  if (msg)   document.getElementById('prog-msg').textContent   = msg;
}

export function showProgress(v) {
  document.getElementById('prog-wrap').style.display = v ? 'block' : 'none';
}

// ── File loaded ──────────────────────────────────────────────
export function setFileLoaded(type, name, stats) {
  document.getElementById(`fn-${type}`).textContent = name;
  document.getElementById(`st-${type}`).textContent = stats;
  document.getElementById(`dz-${type}`).classList.add('ok');
}

// ── Pathway selectors ────────────────────────────────────────
export function populateSelectors(list) {
  const s1 = document.getElementById('sel-single');
  const s2 = document.getElementById('sel-multi');
  s1.innerHTML = '<option value="">— choose pathway —</option>';
  s2.innerHTML = '';
  for (const p of list) {
    const lbl = `${p.name}  (${p.size})`;
    s1.appendChild(new Option(lbl, p.name));
    s2.appendChild(new Option(lbl, p.name));
  }
}

// ── Mode tabs ────────────────────────────────────────────────
export function setupModeTabs() {
  document.querySelectorAll('.pm-tab').forEach(tab => {
    tab.addEventListener('click', () => {
      document.querySelectorAll('.pm-tab')
        .forEach(t => t.classList.remove('active'));
      tab.classList.add('active');
      const m = tab.dataset.mode;
      document.getElementById('sel-single-wrap').style.display =
        m === 'single' ? 'block' : 'none';
      document.getElementById('sel-multi-wrap').style.display  =
        m === 'multi'  ? 'block' : 'none';
    });
  });
}

// ── Get selected pathways ────────────────────────────────────
export function getSelectedPathways(pathwayList) {
  const mode =
    document.querySelector('.pm-tab.active')?.dataset.mode ?? 'all';
  if (mode === 'all') return pathwayList;
  if (mode === 'single') {
    const v = document.getElementById('sel-single').value;
    return v ? pathwayList.filter(p => p.name === v) : pathwayList;
  }
  if (mode === 'multi') {
    const sel  = document.getElementById('sel-multi');
    const vals = new Set(Array.from(sel.selectedOptions, o => o.value));
    return vals.size > 0
      ? pathwayList.filter(p => vals.has(p.name))
      : pathwayList;
  }
  return pathwayList;
}

// ── Run button ───────────────────────────────────────────────
export function setRunEnabled(v) {
  document.getElementById('btn-run').disabled = !v;
}


// ── CSV export ───────────────────────────────────────────────
export function downloadCSV(results, engine) {
  if (!results?.length) return;

  const hdr = [
    '#',           // rank
    'Pathway',     // name
    'Size',        // size
    'NES',         // nes
    'NES_AD',      // nes_ad
    'p_KS',        // pKS
    'p_AD',        // pAD
    'p_Cauchy',    // pCauchy
    'FDR_KS',      // fdr_ks
    'FDR_AD',      // fdr_ad
    'pAD_perm',    // pAD_emp
    'pAD_par'      // pAD_par
  ];

  const fp = v =>
    v == null || isNaN(v) ? '' :
    Math.abs(v) < 0.001 ? v.toExponential(6) : v.toFixed(6);

  const fn = v =>
    v == null || isNaN(v) ? '' : v.toFixed(6);

  const rows = [hdr.join(',')];

  results.forEach((r, idx) => {
    const safeName = `"${r.name.replace(/"/g, '""')}"`;

    const row = [
      idx + 1,             // # (rank)
      safeName,            // Pathway
      r.size ?? '',        // Size
      fn(r.nes),           // NES
      fn(r.nes_ad),        // NES-AD
      fp(r.pKS),           // p_KS
      fp(r.pAD),           // p_AD
      fp(r.pCauchy),       // p_Cauchy
      fp(r.fdr_ks),        // FDR_KS
      fp(r.fdr_ad),        // FDR_AD
      fp(r.pAD_emp),       // pAD (perm)
      fp(r.pAD_par)        // pAD (par)
    ];
    
    rows.push(row.join(','));
  });

  const csvContent = '\uFEFF' + rows.join('\n');
  const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });

  const link = document.createElement('a');
  link.href = URL.createObjectURL(blob);
  link.setAttribute('download', `ingsea_results.csv`);
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}

// ── Demo data ────────────────────────────────────────────────
export function generateDemo() {
  const NG = 200, NS = 20, NC = 10;

  const gNames = Array.from({ length: NG }, (_, i) =>
    `GENE${String(i + 1).padStart(3, '0')}`);
  const sNames = [
    ...Array.from({ length: NC },      (_, i) => `Case_${i + 1}`),
    ...Array.from({ length: NS - NC }, (_, i) => `Ctrl_${i + 1}`)
  ];

  const mat = [];
  for (let g = 0; g < NG; g++) {
    const row  = new Float64Array(NS);
    const base = 4 + Math.random() * 4;
    const sd   = 0.5 + Math.random() * 0.5;
    for (let s = 0; s < NS; s++) {
      let v = base + sd * _randn();
      if (s < NC) {
        if (g < 20)            v += 1.8 + Math.random() * 0.8;
        if (g >= 50 && g < 65) v -= 1.8 + Math.random() * 0.8;
      }
      row[s] = Math.max(0, v);
    }
    mat.push(row);
  }

  const rawPathways = [
    { name: 'DEMO_UPREGULATED',   url: null,
      genes: Array.from({ length: 20 }, (_, i) => gNames[i]) },
    { name: 'DEMO_DOWNREGULATED', url: null,
      genes: Array.from({ length: 15 }, (_, i) => gNames[i + 50]) },
    { name: 'DEMO_MIXED',         url: null,
      genes: [
        ...Array.from({ length: 10 }, (_, i) => gNames[i + 100]),
        ...Array.from({ length: 10 }, (_, i) => gNames[i + 30])
      ]}
  ];

  const pathwayList = buildMasks(rawPathways, gNames);
  return { gNames, sNames, mat, rawPathways, pathwayList, nCase: NC };
}

function _randn() {
  const u = 1 - Math.random();
  const v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
