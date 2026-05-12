// ═══════════════════════════════════════════════════════════
//  main.js  ·  Application entry point
// ═══════════════════════════════════════════════════════════
'use strict';

import { runGSEA }            from './core/gsea.js';
import { parseExpr, parseGMT, parsePathwayCSV, buildMasks, setupDropZone }
  from './ui/fileio.js';
import { drawEnrichmentSVG, updatePlotHeader, renderESStats }
  from './ui/plot.js';
import { renderTable }        from './ui/table.js';
import {
  log, setProgress, showProgress,
  setFileLoaded, populatePathwaySelectors, getSelectedPathways,
  updateRunButton, setupPathwayModeTabs, downloadCSV, generateDemo
} from './ui/controls.js';

// ── Application state ────────────────────────────────────────
const state = {
  exprMat:      null,
  geneNames:    [],
  sampleNames:  [],
  rawPathways:  [],
  pathwayList:  [],
  lastResults:  null,
  engine:       'gg',
  running:      false
};

// ── Initialise UI ────────────────────────────────────────────
setupPathwayModeTabs();

// Extended columns toggle
document.getElementById('chk-extra-cols').addEventListener('change', e => {
  document.getElementById('rt').classList.toggle('show-extra', e.target.checked);
});

// Engine selector
document.getElementById('sel-engine').addEventListener('change', e => {
  state.engine = e.target.value;
});

// ── File drop zones ──────────────────────────────────────────
setupDropZone('dz-expr', 'fi-expr', file => handleExprFile(file));
setupDropZone('dz-path', 'fi-path', file => handlePathFile(file));

function handleExprFile(file) {
  const reader = new FileReader();
  reader.onload = e => {
    try {
      const { gNames, sNames, mat, maxRaw, transformed } = parseExpr(e.target.result);
      state.geneNames   = gNames;
      state.sampleNames = sNames;
      state.exprMat     = mat;

      const note = transformed ? ` [log₂+1, rawMax=${maxRaw.toFixed(0)}]` : '';
      setFileLoaded('expr', file.name, `${gNames.length} genes × ${sNames.length} samples${note}`);
      log(`Expression: ${gNames.length}×${sNames.length}${note}`, 'ok');

      // Update n-case default
      const ncInput = document.getElementById('n-case');
      if (+ncInput.value === 10 && sNames.length !== 20)
        ncInput.value = Math.floor(sNames.length / 2);

      rebuildMasks();
      updateRunButton(!state.running && state.pathwayList.length > 0 && !!state.exprMat);
    } catch (err) {
      log(`Expression error: ${err.message}`, 'err');
    }
  };
  reader.readAsText(file);
}

function handlePathFile(file) {
  const reader = new FileReader();
  reader.onload = e => {
    try {
      const txt = e.target.result;
      let rp;
      const low = file.name.toLowerCase();
      if (low.endsWith('.gmt') || (low.endsWith('.txt') && txt.split('\t').length > 5)) {
        rp = parseGMT(txt);
        if (!rp.length) throw new Error('No valid GMT rows');
      } else {
        rp = parsePathwayCSV(txt);
      }
      state.rawPathways = rp;
      log(`Loaded ${rp.length} gene sets from ${file.name}`, 'ok');

      rebuildMasks();
      const valid = state.pathwayList.length;
      setFileLoaded('path', file.name,
        state.geneNames.length > 0
          ? `${valid} / ${rp.length} sets match expression genes (≥5)`
          : `${rp.length} sets (load expression to match)`
      );
      updateRunButton(!state.running && valid > 0 && !!state.exprMat);
    } catch (err) {
      log(`Gene-set error: ${err.message}`, 'err');
    }
  };
  reader.readAsText(file);
}

function rebuildMasks() {
  if (!state.geneNames.length || !state.rawPathways.length) return;
  state.pathwayList = buildMasks(state.rawPathways, state.geneNames);
  log(`Pathway masks: ${state.pathwayList.length} valid sets (≥5 genes matched)`, 'ok');
  populatePathwaySelectors(state.pathwayList);
}

// ── Run button ───────────────────────────────────────────────
document.getElementById('btn-run').addEventListener('click', async () => {
  if (state.running) return;

  const paths = getSelectedPathways(state.pathwayList);
  if (!paths.length) { log('No pathways selected', 'err'); return; }

  const nCase = +document.getElementById('n-case').value;
  const nPerms = Math.min(+document.getElementById('n-perms').value, 2000);
  const engine = document.getElementById('sel-engine').value;
  state.engine = engine;

  state.running = true;
  updateRunButton(false);
  showProgress(true);
  setProgress(0, 'Starting…');

  log(`Running GSEA: ${paths.length} pathway(s) · ${nPerms} perms · engine=${engine}`, 'ok');

  const t0 = performance.now();
  try {
    const results = await runGSEA({
      exprMat:    state.exprMat,
      geneNames:  state.geneNames,
      pathways:   paths,
      nCase,
      nPerms,
      engine,
      onProgress: (pct, msg) => { setProgress(pct, msg); }
    });

    const elapsed = ((performance.now() - t0) / 1000).toFixed(1);
    log(`Completed in ${elapsed}s · ${results.length} result(s)`, 'ok');
    setProgress(100, `Done in ${elapsed}s`);

    state.lastResults = results;
    renderResults(results, engine);

  } catch (err) {
    log(`Error: ${err.message}`, 'err');
    setProgress(0, `Error: ${err.message}`);
  }

  state.running = false;
  updateRunButton(true);
});

// ── Demo data ────────────────────────────────────────────────
document.getElementById('btn-demo').addEventListener('click', () => {
  const { gNames, sNames, mat, rawPathways, pathwayList, nCase } = generateDemo();
  state.geneNames   = gNames;
  state.sampleNames = sNames;
  state.exprMat     = mat;
  state.rawPathways = rawPathways;
  state.pathwayList = pathwayList;

  document.getElementById('n-case').value = nCase;
  setFileLoaded('expr', 'demo_expression.csv', `${gNames.length} genes × ${sNames.length} samples`);
  setFileLoaded('path', 'demo_pathways', `${pathwayList.length} gene sets`);
  populatePathwaySelectors(pathwayList);
  updateRunButton(true);
  log(`Demo loaded: ${gNames.length}×${sNames.length}, ${pathwayList.length} pathways`, 'ok');
});

// ── Download CSV ─────────────────────────────────────────────
document.getElementById('btn-dl').addEventListener('click', () => {
  downloadCSV(state.lastResults, state.engine);
});

// ── Results renderer ─────────────────────────────────────────
function renderResults(results, engine) {
  document.getElementById('empty-state').style.display  = 'none';
  document.getElementById('res-panel').style.display    = 'block';

  const showFDR = results.length > 10;
  document.getElementById('res-count').textContent =
    `${results.length} pathway${results.length !== 1 ? 's' : ''}`;
  document.getElementById('res-engine').textContent =
    engine === 'gg' ? 'Engine: Γ (KS) + GΓ (AD)' : 'Engine: Empirical only';

  renderTable(results, showFDR, engine, selectPathway);
}

function selectPathway(result) {
  updatePlotHeader(result);
  const container = document.getElementById('svg-container');
  drawEnrichmentSVG(result, state.pathwayList, state.geneNames, container);
  renderESStats(result, state.engine);
}

// ── Resize → redraw SVG ──────────────────────────────────────
let _resizeTimer = null;
window.addEventListener('resize', () => {
  clearTimeout(_resizeTimer);
  _resizeTimer = setTimeout(() => {
    if (!state.lastResults) return;
    const sel = document.querySelector('#tbody tr.sel');
    if (!sel) return;
    const name = sel.dataset.name;
    const r = state.lastResults.find(r => r.name === name);
    if (r) {
      drawEnrichmentSVG(r, state.pathwayList, state.geneNames,
        document.getElementById('svg-container'));
    }
  }, 150);
});

// ── Init log ─────────────────────────────────────────────────
log('iGSEA v2.0 ready — load files or click ⚡ Demo', 'ok');
