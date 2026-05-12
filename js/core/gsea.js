// ═══════════════════════════════════════════════════════════
//  core/gsea.js  ·  GSEA orchestrator
//  All heavy math is in math.js; distribution fitting
//  is delegated to webr-bridge.js.
// ═══════════════════════════════════════════════════════════
'use strict';

import {
  calcSNR, rankOrder, calcGSEAStats,
  cauchyCombine, bhFDR, shuffle,
  calcNES, calcNES_AD, empP_KS, empP_AD
} from './math.js';

import { isWebRReady, fitMany } from './webr-bridge.js';

const MAX_PERMS = 2000;
const CHUNK     = 25;   // permutations per async yield

/**
 * Run GSEA for a set of pathways.
 *
 * @param {object} opts
 *   .exprMat    {Float64Array[]}  genes × samples
 *   .geneNames  {string[]}
 *   .pathways   [{name, mask:Uint8Array, size, url?}]
 *   .nCase      {number}
 *   .nPerms     {number}  clamped to MAX_PERMS
 *   .engine     {'gg'|'empirical'}
 *   .onProgress {(pct:number, phase:string, msg:string) => void}
 * @returns {Promise<Result[]>}
 */
export async function runGSEA(opts) {
  const { exprMat, geneNames, pathways, engine, onProgress } = opts;
  const nPerms = Math.min(opts.nPerms | 0, MAX_PERMS);
  const nG     = geneNames.length;
  const nS     = exprMat[0].length;
  const nCase  = opts.nCase | 0;
  const nP     = pathways.length;

  if (nCase < 2 || nCase > nS - 2)
    throw new Error(
      `nCase=${nCase} is out of range. Need 2 ≤ nCase ≤ ${nS - 2} (nSamples=${nS}).`
    );

  // Fixed case/ctrl index arrays (never change — only permArr is shuffled)
  const cIdx = new Uint16Array(nCase);
  const tIdx = new Uint16Array(nS - nCase);
  for (let i = 0; i < nCase; i++)      cIdx[i] = i;
  for (let i = 0; i < nS - nCase; i++) tIdx[i] = i + nCase;

  // ── 1. Observed statistics ──────────────────────────────────
  onProgress(0, 'Observed', 'Computing observed SNR & enrichment…');
  const obsSnr  = calcSNR(exprMat, nG, cIdx, tIdx);
  const obsOrd  = rankOrder(obsSnr, nG);
  const obsStats = pathways.map(p =>
    calcGSEAStats(obsSnr, obsOrd, p.mask, nG, /*wantCurve=*/true)
  );

  // ── 2. Permutation null distributions ──────────────────────
  // Pre-allocate: one Float64Array per pathway per statistic
  const nullKS = Array.from({ length: nP }, () => new Float64Array(nPerms));
  const nullAD = Array.from({ length: nP }, () => new Float64Array(nPerms));

  // Reusable buffers — allocate ONCE outside the loop
  const permArr = new Uint16Array(nS);
  for (let i = 0; i < nS; i++) permArr[i] = i;
  const permOrd = new Int32Array(nG);

  // Reusable case/ctrl index views into permArr
  // NOTE: subarray returns a VIEW — no copy, always up to date after shuffle
  const permCIdx = permArr.subarray(0, nCase);
  const permTIdx = permArr.subarray(nCase);

  const t0 = performance.now();

  for (let s = 0; s < nPerms; s += CHUNK) {
    const end = Math.min(s + CHUNK, nPerms);

    for (let j = s; j < end; j++) {
      shuffle(permArr);   // modifies permArr in-place → views update automatically

      const pSnr = calcSNR(exprMat, nG, permCIdx, permTIdx);

      // Argsort pSnr descending into reusable permOrd
      for (let i = 0; i < nG; i++) permOrd[i] = i;
      permOrd.sort((a, b) => pSnr[b] - pSnr[a]);

      for (let pi = 0; pi < nP; pi++) {
        const st = calcGSEAStats(pSnr, permOrd, pathways[pi].mask, nG, /*wantCurve=*/false);
        nullKS[pi][j] = st.ks;
        nullAD[pi][j] = st.ad;
      }
    }

    // Progress update + yield
    const elapsed = (performance.now() - t0) / 1000;
    const pct     = Math.round(end / nPerms * 80);   // perms = 0–80%
    const eta     = elapsed > 0 ? ((nPerms - end) / (end / elapsed)).toFixed(0) : '–';
    const rate    = (end / (elapsed + 1e-6)).toFixed(0);
    onProgress(pct, 'Permutations',
      `${end} / ${nPerms} · ETA ${eta}s · ${rate} perm/s`);

    await _yield();
  }

  // ── 3. Empirical p-values & NES (all typed-array, no alloc) ─
  onProgress(81, 'Statistics', 'Computing empirical p-values & NES…');
  await _yield();

  const empResults = pathways.map((p, pi) => {
    const ok   = obsStats[pi].ks;
    const oa   = obsStats[pi].ad;
    const nkA  = nullKS[pi];
    const naA  = nullAD[pi];

    const pKS_emp = empP_KS(ok, nkA, nPerms);
    const pAD_emp = empP_AD(oa, naA, nPerms);
    const nes     = calcNES(ok, nkA);
    const nes_ad  = calcNES_AD(oa, naA);

    return { pi, ok, oa, pKS_emp, pAD_emp, nes, nes_ad };
  });

  // ── 4. Parametric fitting (WebR) ────────────────────────────
  let fitResults;
  const useGG = engine === 'gg' && isWebRReady();

  if (useGG) {
    onProgress(85, 'Fitting', 'Fitting null distributions (R/flexsurv)…');
    await _yield();

    const items = empResults.map((r, pi) => ({
      nullKS: nullKS[pi],
      nullAD: nullAD[pi],
      obsKS:  r.ok,
      obsAD:  r.oa,
      empKS:  r.pKS_emp,
      empAD:  r.pAD_emp
    }));

    fitResults = await fitMany(items);
  } else {
    fitResults = empResults.map(r => ({
      pKS:    r.pKS_emp,
      pAD:    r.pAD_emp,
      engine: engine === 'gg' ? 'empirical-fallback' : 'empirical'
    }));
  }

  // ── 5. Assemble results ─────────────────────────────────────
  onProgress(97, 'Assembling', 'Computing Cauchy combination & FDR…');
  await _yield();

  const results = pathways.map((p, pi) => {
    const er  = empResults[pi];
    const fr  = fitResults[pi];
    const pKS = fr.pKS;
    const pAD = fr.pAD;
    const pC  = cauchyCombine(pKS, pAD);

    return {
      name:     p.name,
      url:      p.url  ?? null,
      size:     p.size,
      es:       er.ok,
      ad:       er.oa,
      nes:      er.nes,
      nes_ad:   er.nes_ad,
      pKS_emp:  er.pKS_emp,
      pAD_emp:  er.pAD_emp,
      pKS:      pKS,
      pAD:      pAD,
      pCauchy:  pC,
      engine:   fr.engine,
      fdr:      null,     // filled below
      curve:    obsStats[pi].curve,
      obsOrd:   obsOrd,   // shared reference — DO NOT MUTATE
      peakIdx:  obsStats[pi].peakIdx
    };
  });

  // BH-FDR on primary p (Cauchy)
  if (results.length > 1) {
    const pvec = results.map(r => r.pCauchy);
    const qvec = bhFDR(pvec);
    results.forEach((r, i) => { r.fdr = qvec[i]; });
  }

  onProgress(100, 'Done', `${results.length} pathway(s) — analysis complete`);
  return results;
}

const _yield = () => new Promise(r => setTimeout(r, 0));
