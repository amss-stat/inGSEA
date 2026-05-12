// ═══════════════════════════════════════════════════════════
//  ui/controls.js  ·  v2.7
//  Changes from v2.6:
//  • CSV export: removed pKS_par column (always null)
//  • CSV header clarifies pKS is always permutation-based
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
  const isP = engine === 'parametric';

  const hdr = [
    'Pathway', 'Size',
    'NES_KS', 'NES_AD',
    'ES', 'AD',
    'pKS_perm',                     // KS is always permutation
    'pAD',                          // combined (par if available, else perm)
    ...(isP ? ['pAD_par'] : []),    // parametric AD only when engine=parametric
    'pAD_perm',                     // always empirical AD
    'pCauchy', 'FDR'
  ];

  const fp = v =>
    v == null ? '' :
    Math.abs(v) < 0.001 ? v.toExponential(6) : v.toFixed(6);

  const rows = [hdr.join(',')];
  for (const r of results) {
    rows.push([
      `"${r.name}"`,
      r.size,
      r.nes.toFixed(6),
      r.nes_ad.toFixed(6),
      r.es.toFixed(6),
      r.ad.toFixed(4),
      fp(r.pKS),                    // always = pKS_emp
      fp(r.pAD),                    // par if fitted, else emp
      ...(isP ? [fp(r.pAD_par)] : []),
      fp(r.pAD_emp),
      fp(r.pCauchy),
      fp(r.fdr)
    ].join(','));
  }

  const blob = new Blob([rows.join('\n')], { type: 'text/csv' });
  const a    = Object.assign(document.createElement('a'), {
    href:     URL.createObjectURL(blob),
    download: 'igsea_results.csv'
  });
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 5000);
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
