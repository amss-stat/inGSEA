// ═══════════════════════════════════════════════════════════
//  main.js  ·  iGSEA application entry
// ═══════════════════════════════════════════════════════════
'use strict';

import { runGSEA }          from './core/gsea.js';
import { initWebR }         from './core/webr-bridge.js';
import {
  parseExpr, parseGMT, buildMasks, setupDropZone
}                           from './ui/fileio.js';
import {
  drawCurve, updatePlotHeader, renderESStats
}                           from './ui/plot.js';
import { renderTable }      from './ui/table.js';
import {
  log, setProgress, showProgress, setFileLoaded,
  populateSelectors, setupModeTabs, getSelectedPathways,
  setRunEnabled, setWebRStatus, downloadCSV, generateDemo
}                           from './ui/controls.js';

// ── State ────────────────────────────────────────────────────
const S = {
  exprMat:     null,
  geneNames:   [],
  sampleNames: [],
  rawPathways: [],
  pathwayList: [],
  results:     null,
  engine:      'parametric',
  running:     false,
  showFDR:     false
};

// ── WebR (background, non-blocking) ──────────────────────────
setWebRStatus('loading');
initWebR(status => {
  setWebRStatus(status);
  if (status === 'ready') log('R engine ready (pure base R, no external packages)', 'ok');
  if (status === 'error') log('R engine unavailable — parametric mode will fall back to permutation', 'warn');
});

// ── UI setup ─────────────────────────────────────────────────
setupModeTabs();

document.getElementById('sel-engine').addEventListener('change', e => {
  S.engine = e.target.value;
});

document.getElementById('chk-extra').addEventListener('change', e => {
  document.getElementById('rt').classList.toggle('show-ext', e.target.checked);
});

// ── File drops ───────────────────────────────────────────────
setupDropZone('dz-expr', 'fi-expr', _loadExpr);
setupDropZone('dz-path', 'fi-path', _loadPath);

function _loadExpr(file) {
  const reader = new FileReader();
  reader.onload = e => {
    try {
      const { gNames, sNames, mat, maxRaw, transformed } = parseExpr(e.target.result);
      S.geneNames   = gNames;
      S.sampleNames = sNames;
      S.exprMat     = mat;

      const note = transformed
        ? ` · log₂+1 (rawMax=${maxRaw.toFixed(0)})` : '';
      setFileLoaded('expr', file.name,
        `${gNames.length} genes × ${sNames.length} samples${note}`);
      log(`Expression: ${gNames.length}×${sNames.length}${note}`, 'ok');

      // Smart nCase default
      const nc = document.getElementById('n-case');
      if (+nc.value === 10 && sNames.length !== 20)
        nc.value = Math.floor(sNames.length / 2);

      _rebuildMasks();
      _updateRun();
    } catch (err) { log(`Expression error: ${err.message}`, 'err'); }
  };
  reader.readAsText(file);
}

function _loadPath(file) {
  const reader = new FileReader();
  reader.onload = e => {
    try {
      // Auto-detect GMT by checking if first line has ≥3 tab-separated fields
      const txt   = e.target.result;
      const first = txt.split('\n')[0] || '';
      const nTabs = (first.match(/\t/g) ?? []).length;

      let rp;
      if (nTabs >= 2) {
        rp = parseGMT(txt);
      } else {
        throw new Error(
          'File does not appear to be GMT format.\n' +
          'Expected tab-separated: name⟨tab⟩url⟨tab⟩gene1⟨tab⟩…\n' +
          'Download GMT files from MSigDB Collections.'
        );
      }

      S.rawPathways = rp;
      log(`Gene sets: ${rp.length} loaded from ${file.name}`, 'ok');
      _rebuildMasks();

      const nV = S.pathwayList.length;
      setFileLoaded('path', file.name,
        S.geneNames.length > 0
          ? `${nV} / ${rp.length} sets (≥10 matched genes)`
          : `${rp.length} sets — load expression to match`
      );
      _updateRun();
    } catch (err) { log(`Gene sets error: ${err.message}`, 'err'); }
  };
  reader.readAsText(file);
}

function _rebuildMasks() {
  if (!S.geneNames.length || !S.rawPathways.length) return;
  S.pathwayList = buildMasks(S.rawPathways, S.geneNames);
  log(`Mask build: ${S.pathwayList.length} pathways with ≥10 matched genes`, 'ok');
  populateSelectors(S.pathwayList);
}

function _updateRun() {
  setRunEnabled(!S.running && !!S.exprMat && S.pathwayList.length > 0);
}

// ── Run iGSEA ────────────────────────────────────────────────
document.getElementById('btn-run').addEventListener('click', async () => {
  if (S.running) return;

  const paths = getSelectedPathways(S.pathwayList);
  if (!paths.length) { log('No pathways selected', 'err'); return; }

  const nCase  = +document.getElementById('n-case').value;
  const nPerms = Math.min(+document.getElementById('n-perms').value, 2000);
  S.engine  = document.getElementById('sel-engine').value;
  S.running = true;
  _updateRun();
  showProgress(true);
  setProgress(0, 'Starting', 'Initialising iGSEA…');

  log(`iGSEA: ${paths.length} pathway(s) · ${nPerms} perms · ` +
      `engine=${S.engine}`, 'ok');
  const t0 = performance.now();

  try {
    const results = await runGSEA({
      exprMat:   S.exprMat,
      geneNames: S.geneNames,
      pathways:  paths,
      nCase,
      nPerms,
      engine:    S.engine,
      onProgress: setProgress
    });

    const sec = ((performance.now() - t0) / 1000).toFixed(1);
    setProgress(100, 'Done', `iGSEA complete in ${sec}s`);
    log(`iGSEA complete: ${results.length} pathway(s) in ${sec}s`, 'ok');

    S.results = results;
    S.showFDR = results.length > 10;
    _renderResults(results);

  } catch (err) {
    log(`Error: ${err.message}`, 'err');
    setProgress(0, 'Error', err.message);
  }

  S.running = false;
  _updateRun();
});

// ── Demo ─────────────────────────────────────────────────────
document.getElementById('btn-demo').addEventListener('click', () => {
  const d = generateDemo();
  Object.assign(S, {
    geneNames:   d.gNames,
    sampleNames: d.sNames,
    exprMat:     d.mat,
    rawPathways: d.rawPathways,
    pathwayList: d.pathwayList
  });
  document.getElementById('n-case').value = d.nCase;
  setFileLoaded('expr', 'demo_expression.csv',
    `${d.gNames.length} genes × ${d.sNames.length} samples`);
  setFileLoaded('path', 'demo_pathways.gmt',
    `${d.pathwayList.length} synthetic gene sets`);
  populateSelectors(d.pathwayList);
  _updateRun();
  log(`Demo loaded: ${d.gNames.length}×${d.sNames.length}, ` +
      `${d.pathwayList.length} pathways`, 'ok');
});

// ── Export ────────────────────────────────────────────────────
document.getElementById('btn-dl').addEventListener('click', () => {
  downloadCSV(S.results, S.engine);
});

// ── Render results ───────────────────────────────────────────
function _renderResults(results) {
  document.getElementById('empty-state').style.display = 'none';
  document.getElementById('res-panel').style.display   = 'block';
  document.getElementById('res-count').textContent =
    `${results.length} pathway${results.length !== 1 ? 's' : ''}`;
  document.getElementById('res-engine-badge').textContent =
    S.engine === 'parametric'
      ? 'Parametric Approximation (Γ + GΓ)'
      : 'Permutation';

  renderTable(results, S.showFDR, S.engine, _selectPathway);
}

function _selectPathway(result) {
  updatePlotHeader(result);
  drawCurve(result, S.pathwayList, S.geneNames,
    document.getElementById('svg-wrap'));
  renderESStats(result, S.engine, S.showFDR);
}

// ── Resize ───────────────────────────────────────────────────
let _rt;
window.addEventListener('resize', () => {
  clearTimeout(_rt);
  _rt = setTimeout(() => {
    if (!S.results) return;
    const sel = document.querySelector('#tbody tr.sel');
    if (!sel) return;
    const r = S.results.find(x => x.name === sel.dataset.name);
    if (r) drawCurve(r, S.pathwayList, S.geneNames,
      document.getElementById('svg-wrap'));
  }, 200);
}, { passive: true });

// ── Init ─────────────────────────────────────────────────────
log('iGSEA v2.2 ready — load files or click ⚡ Demo', 'ok');
