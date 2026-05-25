// ═══════════════════════════════════════════════════════════
//  main.js  ·  adgsea v2.6
// ═══════════════════════════════════════════════════════════
'use strict';

import { runGSEA }          from './core/gsea.js';
import {
  parseExpr, parseGMT, buildMasks, setupDropZone
}                           from './ui/fileio.js';
import {
  drawCurve, updatePlotHeader, renderESStats, resetZoom,
  drawNESChart, exportCurrentCurvePNG, exportAllCurves,
  exportSVGAsPNG
}                           from './ui/plot.js';
import { renderTable }      from './ui/table.js';
import {
  log, setProgress, showProgress, setFileLoaded,
  populateSelectors, setupModeTabs, getSelectedPathways,
  setRunEnabled, downloadCSV, generateDemo
}                           from './ui/controls.js';

const S = {
  exprMat:     null,
  geneNames:   [],
  sampleNames: [],
  rawPathways: [],
  pathwayList: [],
  results:     null,
  engine:      'parametric',
  running:     false,
  showFDR:     false,
  abortCtrl:   null
};

// ── jStat check ──────────────────────────────────────────────
// ── main.js 中的修改 ──────────────────────────────────────────────
function checkJStat() {
  const badge = document.getElementById('jstat-status');
  if (typeof jStat !== 'undefined') {
    log('jStat (Local) ready', 'ok');
    if (badge) {
      badge.textContent = 'jStat ✓';
      badge.className   = 'jstat-badge ok';
    }
  } else {
    log('jStat failed to load locally — check file path', 'err');
    if (badge) {
      badge.textContent = 'jStat ✗';
      badge.className   = 'jstat-badge err';
    }
    const opt = document.querySelector('#sel-engine option[value="parametric"]');
    if (opt) opt.disabled = true;
    document.getElementById('sel-engine').value = 'permutation';
    S.engine = 'permutation';
  }
}
// 脚本一旦执行，本地的 jstat.min.js 肯定已经加载完毕（因为 script 标签阻塞特性）
// 直接同步调用即可：
checkJStat();


// ── Mode tabs ─────────────────────────────────────────────────
setupModeTabs();

// ── Engine selector ───────────────────────────────────────────
document.getElementById('sel-engine').addEventListener('change', e => {
  S.engine = e.target.value;
});

// ── Extended columns toggle ───────────────────────────────────
document.getElementById('chk-extra').addEventListener('change', e => {
  const show = e.target.checked;
  document.getElementById('rt').classList.toggle('show-ext', show);
  // Issue (2): FDR note only visible when extended columns shown AND >10 pathways
  document.getElementById('fdr-note').style.display =
    (show && S.showFDR) ? 'block' : 'none';
});

// ── Reset zoom ────────────────────────────────────────────────
document.getElementById('btn-reset-zoom')
  ?.addEventListener('click', () => resetZoom());

// ── Single curve PNG download ─────────────────────────────────
document.getElementById('btn-dl-one-curve')
  ?.addEventListener('click', () => {
    if (!S.results) return;
    const sel  = document.querySelector('#tbody tr.sel');
    const name = sel?.dataset.name ?? 'curve';
    const safe = name.replace(/[^a-z0-9_\-]/gi, '_').slice(0, 60);
    exportCurrentCurvePNG(`adgsea_${safe}.png`);
  });

// ── Batch curve PNG export ────────────────────────────────────
document.getElementById('btn-dl-curves')
  ?.addEventListener('click', async () => {
    if (!S.results || !S.pathwayList.length) return;
    const btn = document.getElementById('btn-dl-curves');
    btn.disabled    = true;
    btn.textContent = '⏳ Exporting…';
    try {
      await exportAllCurves(S.results, S.pathwayList, S.geneNames);
    } finally {
      btn.disabled    = false;
      btn.textContent = '🖼 Export curves';
    }
  });

// ── NES chart PNG download ────────────────────────────────────
document.getElementById('btn-dl-nes')
  ?.addEventListener('click', () => {
    const svgEl = document.querySelector('#nes-chart-wrap #nes-svg');
    if (!svgEl) return;
    const date = new Date().toISOString().slice(0, 10);
    exportSVGAsPNG(svgEl, `adgsea_nes_chart_${date}.png`);
  });

// ── NES chart toggle (>20 pathways) ──────────────────────────
document.getElementById('btn-toggle-nes')
  ?.addEventListener('click', () => {
    const body = document.getElementById('nes-chart-body');
    const btn  = document.getElementById('btn-toggle-nes');
    const collapsed = body.style.display === 'none';
    body.style.display   = collapsed ? 'block' : 'none';
    btn.textContent      = collapsed ? '▲ Hide chart' : '▼ Show chart';
  });

// ── Clear ─────────────────────────────────────────────────────
document.getElementById('btn-clear').addEventListener('click', () => {
  if (S.running) return;
  S.results = null;
  S.showFDR = false;
  document.getElementById('res-panel').style.display         = 'none';
  document.getElementById('empty-state').style.display       = 'flex';
  document.getElementById('plot-section').style.display      = 'none';
  document.getElementById('nes-chart-section').style.display = 'none';
  document.getElementById('prog-wrap').style.display         = 'none';
  document.getElementById('fdr-note').style.display          = 'none';
  document.getElementById('tbody').innerHTML = '';
  // Reset extended columns checkbox
  const chk = document.getElementById('chk-extra');
  if (chk.checked) {
    chk.checked = false;
    document.getElementById('rt').classList.remove('show-ext');
  }
  log('Results cleared', 'ok');
});

// ── Abort ─────────────────────────────────────────────────────
document.getElementById('btn-abort').addEventListener('click', () => {
  if (!S.running || !S.abortCtrl) return;
  S.abortCtrl.abort();
  document.getElementById('btn-abort').disabled = true;
  log('Abort requested…', 'warn');
});

// ── File drops ────────────────────────────────────────────────
setupDropZone('dz-expr', 'fi-expr', _loadExpr);
setupDropZone('dz-path', 'fi-path', _loadPath);

function _loadExpr(file) {
  const r = new FileReader();
  r.onload = e => {
    try {
      const { gNames, sNames, mat, maxRaw, transformed } =
        parseExpr(e.target.result);
      S.geneNames   = gNames;
      S.sampleNames = sNames;
      S.exprMat     = mat;
      const note = transformed
        ? ` · log₂+1 (rawMax=${maxRaw.toFixed(0)})` : '';
      setFileLoaded('expr', file.name,
        `${gNames.length} genes × ${sNames.length} samples${note}`);
      log(`Expression: ${gNames.length}×${sNames.length}${note}`, 'ok');
      const nc = document.getElementById('n-case');
      if (+nc.value === 10 && sNames.length !== 20)
        nc.value = Math.floor(sNames.length / 2);
      _rebuildMasks();
      _updateRun();
    } catch (err) { log(`Expression error: ${err.message}`, 'err'); }
  };
  r.readAsText(file);
}

function _loadPath(file) {
  const r = new FileReader();
  r.onload = e => {
    try {
      const txt   = e.target.result;
      const first = txt.split('\n')[0] || '';
      if ((first.match(/\t/g) ?? []).length < 2)
        throw new Error(
          'Not GMT format. Expected: name⟨tab⟩url⟨tab⟩gene1⟨tab⟩…'
        );
      const rp = parseGMT(txt);
      S.rawPathways = rp;
      log(`Gene sets: ${rp.length} loaded`, 'ok');
      _rebuildMasks();
      setFileLoaded('path', file.name,
        S.geneNames.length > 0
          ? `${S.pathwayList.length} / ${rp.length} sets (≥10 matched genes)`
          : `${rp.length} sets — load expression to match`);
      _updateRun();
    } catch (err) { log(`Gene sets: ${err.message}`, 'err'); }
  };
  r.readAsText(file);
}

function _rebuildMasks() {
  if (!S.geneNames.length || !S.rawPathways.length) return;
  S.pathwayList = buildMasks(S.rawPathways, S.geneNames);
  log(`Masks: ${S.pathwayList.length} pathways with ≥10 matched genes`, 'ok');
  populateSelectors(S.pathwayList);
}

function _updateRun() {
  setRunEnabled(!S.running && !!S.exprMat && S.pathwayList.length > 0);
  document.getElementById('btn-abort').disabled = !S.running;
  document.getElementById('btn-clear').disabled =  S.running;
}

// ── Run adgsea ─────────────────────────────────────────────────
document.getElementById('btn-run').addEventListener('click', async () => {
  if (S.running) return;

  const paths  = getSelectedPathways(S.pathwayList);
  if (!paths.length) { log('No pathways selected', 'err'); return; }

  const nCase  = +document.getElementById('n-case').value;
  const nPerms = Math.min(+document.getElementById('n-perms').value, 2000);
  const wt     = Math.max(0, +document.getElementById('weight-p').value);
  S.engine    = document.getElementById('sel-engine').value;
  S.running   = true;
  S.abortCtrl = new AbortController();

  _updateRun();
  showProgress(true);
  setProgress(0, 'Starting', 'Initialising adgsea…');
  log(
    `adgsea: ${paths.length} pathway(s) · ${nPerms} perms · ` +
    `p=${wt} · engine=${S.engine}`, 'ok'
  );

  const t0 = performance.now();
  try {
    const results = await runGSEA({
      exprMat:     S.exprMat,
      geneNames:   S.geneNames,
      pathways:    paths,
      nCase,
      nPerms,
      weightP:     wt,
      engine:      S.engine,
      abortSignal: S.abortCtrl.signal,
      onProgress:  setProgress
    });

    const sec = ((performance.now() - t0) / 1000).toFixed(1);
    setProgress(100, 'Done', `adgsea complete in ${sec}s`);
    log(`Done: ${results.length} pathway(s) in ${sec}s`, 'ok');

    S.results = results;
    S.showFDR = results.length > 10;
    _renderResults(results);

  } catch (err) {
    if (err.message === 'Aborted') {
      log('Aborted by user', 'warn');
      setProgress(0, 'Aborted', 'Stopped. Press Run adgsea to restart.');
    } else {
      log(`Error: ${err.message}`, 'err');
      setProgress(0, 'Error', err.message);
    }
  }

  S.running   = false;
  S.abortCtrl = null;
  _updateRun();
});

// ── Demo ──────────────────────────────────────────────────────
document.getElementById('btn-demo').addEventListener('click', () => {
  if (S.running) return;
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
  log(
    `Demo: ${d.gNames.length}×${d.sNames.length}, ` +
    `${d.pathwayList.length} pathways`, 'ok'
  );
});

// ── CSV export ────────────────────────────────────────────────
document.getElementById('btn-dl').addEventListener('click', () => {
  downloadCSV(S.results, S.engine);
});

// ── Render results ────────────────────────────────────────────
function _renderResults(results) {
  document.getElementById('empty-state').style.display = 'none';
  document.getElementById('res-panel').style.display   = 'block';

  // Issue (5): pathway count only
  document.getElementById('res-count').textContent =
    `${results.length} pathway${results.length !== 1 ? 's' : ''}`;

  // Issue (2): FDR note hidden by default; appears only when
  // extended columns are toggled on AND >10 pathways
  document.getElementById('fdr-note').style.display = 'none';

  // Issue (3): NES chart — shown when >10 pathways
  const nesSection   = document.getElementById('nes-chart-section');
  const nesBody      = document.getElementById('nes-chart-body');
  const nesToggleBtn = document.getElementById('btn-toggle-nes');

  if (results.length > 10) {
    nesSection.style.display = 'block';
    drawNESChart(results, document.getElementById('nes-chart-wrap'));

    // Issue (3): collapse by default when >20 pathways
    if (results.length > 20) {
      nesBody.style.display      = 'none';
      nesToggleBtn.style.display = 'inline-flex';
      nesToggleBtn.textContent   = '▼ Show chart';
    } else {
      nesBody.style.display      = 'block';
      nesToggleBtn.style.display = 'none';
    }
  } else {
    nesSection.style.display = 'none';
  }

  renderTable(results, S.showFDR, S.engine, _selectPathway);
}

function _selectPathway(result) {
  updatePlotHeader(result);
  drawCurve(
    result, S.pathwayList, S.geneNames,
    document.getElementById('svg-wrap')
  );
  renderESStats(result, S.engine, S.showFDR);
}

// ── Resize ────────────────────────────────────────────────────
let _resizeTimer;
window.addEventListener('resize', () => {
  clearTimeout(_resizeTimer);
  _resizeTimer = setTimeout(() => {
    // Redraw NES chart if visible
    if (S.results && S.results.length > 10) {
      const wrap = document.getElementById('nes-chart-wrap');
      if (wrap && document.getElementById('nes-chart-body').style.display !== 'none')
        drawNESChart(S.results, wrap);
    }
    // Redraw selected enrichment curve
    if (!S.results) return;
    const sel = document.querySelector('#tbody tr.sel');
    if (!sel) return;
    const r = S.results.find(x => x.name === sel.dataset.name);
    if (r) drawCurve(
      r, S.pathwayList, S.geneNames,
      document.getElementById('svg-wrap')
    );
  }, 200);
}, { passive: true });

log('adgsea v2.6 ready — load files or click ⚡ Demo', 'ok');
