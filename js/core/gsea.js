// ═══════════════════════════════════════════════════════════
//  core/gsea.js  ·  GSEA orchestrator
// ═══════════════════════════════════════════════════════════
'use strict';

import {
  calcSNR, rankOrder, calcGSEAStats,
  cauchyCombine, bhFDR, shuffle,
  calcNES, calcNES_AD, empP_KS, empP_AD
} from './math.js';

import { isWebRReady, fitMany } from './webr-bridge.js';

const MAX_PERMS = 2000;

/**
 * Adaptive chunk size: larger gene sets need smaller chunks
 * to avoid blocking, but small sets can use big chunks for speed.
 * Target: ~40ms per chunk on a mid-range machine.
 */
function chooseChunk(nGenes, nPathways) {
  const cost = nGenes * nPathways;
  if (cost < 20000)   return 100;    // tiny (demo data)
  if (cost < 100000)  return 50;
  if (cost < 500000)  return 25;
  return 10;
}

/**
 * Run iGSEA for a set of pathways.
 */
export async function runGSEA(opts) {
  const { exprMat, geneNames, pathways, engine, onProgress } = opts;
  const nPerms = Math.min(opts.nPerms | 0, MAX_PERMS);
  const nG     = geneNames.length;
  const nS     = exprMat[0].length;
  const nCase  = opts.nCase | 0;
  const nP     = pathways.length;

  if (nCase < 2 || nCase > nS - 2) {
    throw new Error(
      `nCase=${nCase} out of range [2, ${nS - 2}] for ${nS} samples.`
    );
  }

  const CHUNK = chooseChunk(nG, nP);

  // Fixed case/ctrl indices
  const cIdx = new Uint16Array(nCase);
  const tIdx = new Uint16Array(nS - nCase);
  for (let i = 0; i < nCase; i++)      cIdx[i] = i;
  for (let i = 0; i < nS - nCase; i++) tIdx[i] = i + nCase;

  // ── 1. Observed ─────────────────────────────────────────────
  onProgress(0, 'Observed', 'Computing SNR & enrichment statistics…');
  const obsSnr = calcSNR(exprMat, nG, cIdx, tIdx);
  const obsOrd = rankOrder(obsSnr, nG);

  // Observed stats with curves
  const obsStats = pathways.map(p =>
    calcGSEAStats(obsSnr, obsOrd, p.mask, nG, true)
  );

  // ── 2. Permutations ─────────────────────────────────────────
  const nullKS = Array.from({ length: nP }, () => new Float64Array(nPerms));
  const nullAD = Array.from({ length: nP }, () => new Float64Array(nPerms));

  // Reusable buffers
  const permArr  = new Uint16Array(nS);
  for (let i = 0; i < nS; i++) permArr[i] = i;
  const permOrd  = new Int32Array(nG);
  // Views into shuffled array — zero-copy
  const permC = permArr.subarray(0, nCase);
  const permT = permArr.subarray(nCase);

  const t0 = performance.now();

  for (let s = 0; s < nPerms; s += CHUNK) {
    const end = Math.min(s + CHUNK, nPerms);

    for (let j = s; j < end; j++) {
      shuffle(permArr);

      const pSnr = calcSNR(exprMat, nG, permC, permT);

      for (let i = 0; i < nG; i++) permOrd[i] = i;
      permOrd.sort((a, b) => pSnr[b] - pSnr[a]);

      for (let pi = 0; pi < nP; pi++) {
        const st = calcGSEAStats(pSnr, permOrd, pathways[pi].mask, nG, false);
        nullKS[pi][j] = st.ks;
        nullAD[pi][j] = st.ad;
      }
    }

    const elapsed = (performance.now() - t0) / 1000;
    const pct     = Math.round(end / nPerms * 80);
    const eta     = elapsed > 0.1
      ? ((nPerms - end) / (end / elapsed)).toFixed(0) : '—';
    const rate    = (end / (elapsed + 0.001)).toFixed(0);
    onProgress(pct, 'Permutations',
      `${end}/${nPerms} · ETA ${eta}s · ${rate} perm/s`);

    await _yield();
  }

  // ── 3. Empirical p & NES ────────────────────────────────────
  onProgress(81, 'Statistics', 'Computing empirical p-values & NES…');
  await _yield();

  const empArr = [];
  for (let pi = 0; pi < nP; pi++) {
    const ok = obsStats[pi].ks;
    const oa = obsStats[pi].ad;
    empArr.push({
      ok, oa,
      pKS_emp: empP_KS(ok, nullKS[pi], nPerms),
      pAD_emp: empP_AD(oa, nullAD[pi], nPerms),
      nes:     calcNES(ok, nullKS[pi]),
      nes_ad:  calcNES_AD(oa, nullAD[pi])
    });
  }

  // ── 4. Parametric fitting ───────────────────────────────────
  let fitArr;
  const useParam = engine === 'parametric' && isWebRReady();

  if (useParam) {
    onProgress(85, 'Fitting', 'Fitting null distributions via R…');
    await _yield();

    const items = empArr.map((r, pi) => ({
      nullKS: nullKS[pi],
      nullAD: nullAD[pi],
      obsKS:  r.ok,
      obsAD:  r.oa,
      empKS:  r.pKS_emp,
      empAD:  r.pAD_emp
    }));

    fitArr = await fitMany(items);
  } else {
    fitArr = empArr.map(r => ({
      pKS:    r.pKS_emp,
      pAD:    r.pAD_emp,
      engine: engine === 'parametric' ? 'permutation-fallback' : 'permutation'
    }));
  }

  // ── 5. Assemble results ─────────────────────────────────────
  onProgress(97, 'Assembling', 'Computing Cauchy p & FDR…');
  await _yield();

  const results = pathways.map((p, pi) => {
    const e  = empArr[pi];
    const f  = fitArr[pi];
    const pC = cauchyCombine(f.pKS, f.pAD);

    return {
      name:    p.name,
      url:     p.url ?? null,
      size:    p.size,
      es:      e.ok,
      ad:      e.oa,
      nes:     e.nes,
      nes_ad:  e.nes_ad,
      pKS_emp: e.pKS_emp,
      pAD_emp: e.pAD_emp,
      pKS:     f.pKS,
      pAD:     f.pAD,
      pCauchy: pC,
      engine:  f.engine,
      fdr:     null,
      curve:   obsStats[pi].curve,
      obsOrd:  obsOrd,   // shared
      peakIdx: obsStats[pi].peakIdx
    };
  });

  // BH-FDR
  if (results.length > 1) {
    const pv = results.map(r => r.pCauchy);
    const qv = bhFDR(pv);
    results.forEach((r, i) => { r.fdr = qv[i]; });
  }

  onProgress(100, 'Done', `${results.length} pathway(s) complete`);
  return results;
}

const _yield = () => new Promise(r => setTimeout(r, 0));
