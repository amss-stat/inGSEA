// ═══════════════════════════════════════════════════════════
//  ui/controls.js  ·  UI state helpers & demo-data generator
// ═══════════════════════════════════════════════════════════
'use strict';

import { buildMasks } from './fileio.js';

// ── Logging ─────────────────────────────────────────────────
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
export function setProgress(pct, msg) {
  document.getElementById('prog-fill').style.width = pct + '%';
  document.getElementById('prog-pct').textContent  = pct + '%';
  if (msg) document.getElementById('prog-status').textContent = msg;
}

export function showProgress(visible) {
  document.getElementById('prog-wrap').style.display = visible ? 'block' : 'none';
}

// ── File status ──────────────────────────────────────────────
export function setFileLoaded(type, filename, stats) {
  document.getElementById(`fn-${type}`).textContent  = filename;
  document.getElementById(`st-${type}`).textContent  = stats;
  document.getElementById(`dz-${type}`).classList.add('ok');
}

// ── Pathway selector ─────────────────────────────────────────
export function populatePathwaySelectors(pathwayList) {
  const singleSel = document.getElementById('sel-path-single');
  const multiSel  = document.getElementById('sel-path-multi');

  singleSel.innerHTML = '<option value="">— choose pathway —</option>';
  multiSel.innerHTML  = '';

  for (const p of pathwayList) {
    const label = `${p.name}  (${p.size})`;
    const o1 = new Option(label, p.name);
    const o2 = new Option(label, p.name);
    singleSel.appendChild(o1);
    multiSel.appendChild(o2);
  }
}

/** Return selected pathways based on current mode tab. */
export function getSelectedPathways(pathwayList) {
  const activeTab = document.querySelector('.pm-tab.active');
  const mode = activeTab?.dataset.mode || 'all';

  if (mode === 'all') return pathwayList;

  if (mode === 'single') {
    const val = document.getElementById('sel-path-single').value;
    return val ? pathwayList.filter(p => p.name === val) : pathwayList;
  }

  if (mode === 'multi') {
    const sel = document.getElementById('sel-path-multi');
    const vals = new Set(Array.from(sel.selectedOptions).map(o => o.value));
    return vals.size > 0 ? pathwayList.filter(p => vals.has(p.name)) : pathwayList;
  }

  return pathwayList;
}

// ── Run button ───────────────────────────────────────────────
export function updateRunButton(enabled) {
  document.getElementById('btn-run').disabled = !enabled;
}

// ── Pathway mode tabs ────────────────────────────────────────
export function setupPathwayModeTabs() {
  document.querySelectorAll('.pm-tab').forEach(tab => {
    tab.addEventListener('click', () => {
      document.querySelectorAll('.pm-tab').forEach(t => t.classList.remove('active'));
      tab.classList.add('active');
      const mode = tab.dataset.mode;
      document.getElementById('sel-single-wrap').style.display = mode === 'single' ? 'block' : 'none';
      document.getElementById('sel-multi-wrap').style.display  = mode === 'multi'  ? 'block' : 'none';
    });
  });
}

// ── Download CSV ─────────────────────────────────────────────
export function downloadCSV(results, engine) {
  if (!results?.length) return;
  const rows = [['Pathway','Size','NES_KS','NES_AD','ES','AD',
    engine === 'gg' ? 'pKS_gamma' : 'pKS_emp',
    engine === 'gg' ? 'pAD_GG'    : 'pAD_emp',
    'pCauchy','FDR','pKS_emp','pAD_emp']];
  for (const r of results) {
    rows.push([
      r.name, r.size,
      r.nes.toFixed(4), r.nes_ad.toFixed(4),
      r.es.toFixed(4),  r.ad.toFixed(4),
      r.pKS_fit, r.pAD_fit, r.pCauchy,
      r.fdr ?? '', r.pKS_emp, r.pAD_emp
    ]);
  }
  const csv = rows.map(r => r.join(',')).join('\n');
  const a   = document.createElement('a');
  a.href    = 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv);
  a.download = 'igsea_results.csv';
  a.click();
}

// ── Demo data generator ──────────────────────────────────────
export function generateDemo() {
  const NG = 200, NS = 20, NC = 10;

  const gNames = Array.from({ length: NG }, (_, i) =>
    `GENE${String(i + 1).padStart(3, '0')}`);
  const sNames = [
    ...Array.from({ length: NC },      (_, i) => `Case_${i + 1}`),
    ...Array.from({ length: NS - NC }, (_, i) => `Ctrl_${i + 1}`)
  ];

  const up  = new Set(Array.from({ length: 20 }, (_, i) => i));
  const dn  = new Set(Array.from({ length: 15 }, (_, i) => i + 50));

  const mat = [];
  for (let g = 0; g < NG; g++) {
    const row = new Float64Array(NS);
    const base = 4 + Math.random() * 4;
    const sd   = 0.5 + Math.random() * 0.5;
    for (let s = 0; s < NS; s++) {
      let v = base + sd * randn();
      if (s < NC) {
        if (up.has(g)) v += 1.8 + Math.random() * 0.8;
        if (dn.has(g)) v -= 1.8 + Math.random() * 0.8;
      }
      row[s] = Math.max(0, v);
    }
    mat.push(row);
  }

  const rawPathways = [
    { name: 'DEMO_UPREGULATED',   url: null, genes: Array.from({ length: 20 }, (_, i) => gNames[i]) },
    { name: 'DEMO_DOWNREGULATED', url: null, genes: Array.from({ length: 15 }, (_, i) => gNames[i + 50]) },
    { name: 'DEMO_MIXED',         url: null, genes: [
      ...Array.from({ length: 10 }, (_, i) => gNames[i + 100]),
      ...Array.from({ length: 10 }, (_, i) => gNames[i + 30])
    ]}
  ];

  const pathwayList = buildMasks(rawPathways, gNames);

  return { gNames, sNames, mat, rawPathways, pathwayList, nCase: NC };
}

function randn() {
  const u = 1 - Math.random(), v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
