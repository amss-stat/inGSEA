// ═══════════════════════════════════════════════════════════
//  main.js  ·  Application entry point
// ═══════════════════════════════════════════════════════════
'use strict';

import { runGSEA }          from './core/gsea.js';
import { initWebR }         from './core/webr-bridge.js';
import {
  parseExpr, parseGMT, parsePathwayCSV,
  buildMasks, setupDropZone
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

// ── Application state ────────────────────────────────────────
const S = {
  exprMat:     null,
  geneNames:   [],
  sampleNames: [],
  rawPathways: [],
  pathwayList: [],
  results:     null,
  engine:      'gg',
  running:     false,
  showFDR:     false
};

// ── WebR init (background, non-blocking) ─────────────────────
setWebRStatus('loading');
initWebR(status => {
  setWebRStatus(status);
  if (status === 'ready')  log('WebR + flexsurv ready', 'ok');
  if (status === 'error')  log('WebR failed — will use empirical p-values', 'warn');
  if (status === 'installing') log('Installing R packages…');
});

// ── Mode tabs ────────────────────────────────────────────────
setupModeTabs();

// ── Engine selector ───────────────────────────────────────────
document.getElementById('sel-engine').addEventListener('change', e => {
  S.engine = e.target.value;
});

// ── Extended columns toggle ───────────────────────────────────
document.getElementById('chk-extra').addEventListener('change', e => {
  document.getElementById('rt').classList.toggle('show-ext', e.target.checked);
});

// ── File drop zones ──────────────────────────────────────────
setupDropZone('dz-expr', 'fi-expr', file => _loadExpr(file));
setupDropZone('dz-path', 'fi-path', file => _loadPath(file));

function _loadExpr(file) {
  const r = new FileReader();
  r.onload = e => {
    try {
      const { gNames, sNames, mat, maxRaw, transformed } = parseExpr(e.target.result);
      S.geneNames   = gNames;
      S.sampleNames = sNames;
      S.exprMat     = mat;

      const note = transformed
        ? ` · log₂+1 applied (rawMax=${maxRaw.toFixed(0)})`
        : '';
      setFileLoaded('expr', file.name,
        `${gNames.length} genes × ${sNames.length} samples${note}`);
      log(`Expression: ${gNames.length}×${sNames.length}${note}`, 'ok');

      // Suggest a sensible nCase default
      const nc = document.getElementById('n-case');
      if (+nc.value === 10) nc.value = Math.floor(sNames.length / 2);

      _rebuildMasks();
      _updateRunBtn();
    } catch (err) { log(`Expression: ${err.message}`, 'err'); }
  };
  r.readAsText(file);
}

function _loadPath(file) {
  const r = new FileReader();
  r.onload = e => {
    try {
      const txt = e.target.result;
      const low = file.name.toLowerCase();
      let rp;
      // Detect GMT: tab count per line is high, or extension is .gmt
      if (low.endsWith('.gmt') || txt.split('\n')[0].split('\t').length > 4) {
        rp = parseGMT(txt);
      } else {
        rp = parsePathwayCSV(txt);
      }
      S.rawPathways = rp;
      log(`Gene sets: ${rp.length} loaded from ${file.name}`, 'ok');

      _rebuildMasks();
      const nValid = S.pathwayList.length;
      setFileLoaded('path', file.name,
        S.geneNames.length > 0
          ? `${nValid} / ${rp.length} sets (≥5 matched genes)`
          : `${rp.length} sets (load expression to match)`
      );
      _updateRunBtn();
    } catch (err) { log(`Gene sets: ${err.message}`, 'err'); }
  };
  r.readAsText(file);
}

function _rebuildMasks() {
  if (!S.geneNames.length || !S.rawPathways.length) return;
  S.pathwayList = buildMasks(S.rawPathways, S.geneNames);
  log(`Masks built: ${S.pathwayList.length} valid pathways (≥5 genes matched)`, 'ok');
  populateSelectors(S.pathwayList);
}

function _updateRunBtn() {
  setRunEnabled(!S.running && !!S.exprMat && S.pathwayList.length > 0);
}

// ── Run ──────────────────────────────────────────────────────
document.getElementById('btn-run').addEventListener('click', async () => {
  if (S.running) return;

  const paths = getSelectedPathways(S.pathwayList);
  if (!paths.length) { log('No pathways selected', 'err'); return; }

  const nCase  = +document.getElementById('n-case').value;
  const nPerms = Math.min(+document.getElementById('n-perms').value, 2000);
  S.engine  = document.getElementById('sel-engine').value;
  S.running = true;
  _updateRunBtn();
  showProgress(true);
  setProgress(0, 'Starting', 'Initialising…');

  log(`GSEA start: ${paths.length} pathway(s) · ${nPerms} perms · engine=${S.engine}`, 'ok');
  const t0 = performance.now();

  try {
    const results = await runGSEA({
      exprMat:    S.exprMat,
      geneNames:  S.geneNames,
      pathways:   paths,
      nCase,
      nPerms,
      engine:     S.engine,
      onProgress: (pct, phase, msg) => setProgress(pct, phase, msg)
    });

    const sec = ((performance.now() - t0) / 1000).toFixed(1);
    setProgress(100, 'Done', `${results.length} pathway(s) complete in ${sec}s`);
    log(`Done in ${sec}s`, 'ok');

    S.results  = results;
    S.showFDR  = results.length > 10;
    _renderResults(results);

  } catch (err) {
    log(`Error: ${err.message}`, 'err');
    setProgress(0, 'Error', err.message);
  }

  S.running = false;
  _updateRunBtn();
});

// ── Demo ─────────────────────────────────────────────────────
document.getElementById('btn-demo').addEventListener('click', () => {
  const { gNames, sNames, mat, rawPathways, pathwayList, nCase } = generateDemo();
  Object.assign(S, {
    geneNames:   gNames,
    sampleNames: sNames,
    exprMat:     mat,
    rawPathways,
    pathwayList
  });
  document.getElementById('n-case').value = nCase;
  setFileLoaded('expr', 'demo_expression.csv',
    `${gNames.length} genes × ${sNames.length} samples`);
  setFileLoaded('path', 'demo_pathways.gmt',
    `${pathwayList.length} synthetic gene sets`);
  populateSelectors(pathwayList);
  _updateRunBtn();
  log(`Demo loaded: ${gNames.length}×${sNames.length}, ${pathwayList.length} pathways`, 'ok');
});

// ── CSV export ────────────────────────────────────────────────
document.getElementById('btn-dl').addEventListener('click', () => {
  downloadCSV(S.results, S.engine);
});

// ── Results rendering ─────────────────────────────────────────
function _renderResults(results) {
  document.getElementById('empty-state').style.display = 'none';
  document.getElementById('res-panel').style.display   = 'block';
  document.getElementById('res-count').textContent     =
    `${results.length} pathway${results.length !== 1 ? 's' : ''}`;
  document.getElementById('res-engine-badge').textContent =
    S.engine === 'gg' ? 'Γ (KS) + GΓ (AD) via WebR' : 'Empirical permutation';

  renderTable(results, S.showFDR, S.engine, _selectPathway);
}

function _selectPathway(result) {
  updatePlotHeader(result);
  drawCurve(result, S.pathwayList, S.geneNames,
    document.getElementById('svg-wrap'));
  renderESStats(result, S.engine, S.showFDR);
}

// ── Resize → redraw ───────────────────────────────────────────
let _rTimer = null;
window.addEventListener('resize', () => {
  clearTimeout(_rTimer);
  _rTimer = setTimeout(() => {
    if (!S.results) return;
    const sel  = document.querySelector('#tbody tr.sel');
    if (!sel) return;
    const name = sel.dataset.name;
    const r    = S.results.find(x => x.name === name);
    if (r) drawCurve(r, S.pathwayList, S.geneNames,
      document.getElementById('svg-wrap'));
  }, 150);
}, { passive: true });

// ── Init ─────────────────────────────────────────────────────
log('iGSEA v2.1 ready — load files or click ⚡ Demo', 'ok');
