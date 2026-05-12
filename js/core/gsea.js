// ═══════════════════════════════════════════════════════════
//  core/gsea.js  ·  Async GSEA engine
// ═══════════════════════════════════════════════════════════
'use strict';

import {
  calcSNR, rankOrder, calcGSEAStats,
  cauchyCombine, bhFDR, shuffle, calcNES, calcNES_AD
} from './math.js';

import {
  fitGamma, fitGenGamma, pvalFromGamma, pvalFromGG
} from './distributions.js';

const MAX_PERMS = 2000;
const CHUNK     = 20;   // permutations per async chunk

/**
 * Main GSEA runner.
 *
 * @param {object} opts
 *   .exprMat      Float64Array[]   – genes × samples
 *   .geneNames    string[]
 *   .pathways     [{name,mask,size,url?}]
 *   .nCase        number           – number of case columns
 *   .nPerms       number           – clamped to MAX_PERMS
 *   .engine       'gg'|'empirical'
 *   .onProgress   (pct, msg) => void
 * @returns Promise<Result[]>
 */
export async function runGSEA(opts) {
  const { exprMat, geneNames, pathways, engine, onProgress } = opts;
  const nPerms = Math.min(opts.nPerms, MAX_PERMS);
  const nG = geneNames.length;
  const nS = exprMat[0].length;
  const nCase = opts.nCase;

  if (nCase < 2 || nCase > nS - 2)
    throw new Error(`nCase must be between 2 and ${nS - 2} (total samples: ${nS})`);

  const cIdx = Array.from({ length: nCase },       (_, i) => i);
  const tIdx = Array.from({ length: nS - nCase },  (_, i) => i + nCase);
  const nP   = pathways.length;

  // ── 1. Observed statistics ──────────────────────────────
  onProgress(0, 'Computing observed statistics…');
  const obsSnr = calcSNR(exprMat, nG, cIdx, tIdx);
  const obsOrd = rankOrder(obsSnr, nG);

  const obsStats = pathways.map(p =>
    calcGSEAStats(obsSnr, obsOrd, p.mask, nG, true)
  );

  // ── 2. Permutation null distributions ──────────────────
  const nullKS = Array.from({ length: nP }, () => new Float64Array(nPerms));
  const nullAD = Array.from({ length: nP }, () => new Float64Array(nPerms));

  const permArr = Uint16Array.from({ length: nS }, (_, i) => i);
  const permOrd = new Int32Array(nG);

  const t0 = performance.now();

  for (let s = 0; s < nPerms; s += CHUNK) {
    const end = Math.min(s + CHUNK, nPerms);

    for (let j = s; j < end; j++) {
      shuffle(permArr);

      // Views of permArr for case/ctrl (avoids copy)
      const pSnr = calcSNR(exprMat, nG,
        permArr.subarray(0, nCase),
        permArr.subarray(nCase)
      );
      for (let i = 0; i < nG; i++) permOrd[i] = i;
      permOrd.sort((a, b) => pSnr[b] - pSnr[a]);

      for (let pi = 0; pi < nP; pi++) {
        const st = calcGSEAStats(pSnr, permOrd, pathways[pi].mask, nG, false);
        nullKS[pi][j] = st.ks;
        nullAD[pi][j] = st.ad;
      }
    }

    const elapsed = (performance.now() - t0) / 1000;
    const pct = Math.round(end / nPerms * 100);
    const eta = end > 0 ? ((nPerms - end) / end * elapsed).toFixed(0) : '–';
    const rate = Math.round(end / (elapsed + 1e-6));
    onProgress(pct,
      `Permutation ${end} / ${nPerms} · ETA ${eta}s · ${rate} perm/s`
    );

    await yieldToBrowser();
  }

  // ── 3. Fit parametric distributions & compute p-values ─
  onProgress(98, 'Fitting null distributions…');
  await yieldToBrowser();

  const results = [];

  for (let pi = 0; pi < nP; pi++) {
    const p      = pathways[pi];
    const ok     = obsStats[pi].ks;
    const oa     = obsStats[pi].ad;
    const nkArr  = nullKS[pi];
    const naArr  = nullAD[pi];

    // Empirical p-values
    let ck = 0, ca = 0;
    for (let j = 0; j < nPerms; j++) {
      if (Math.abs(nkArr[j]) >= Math.abs(ok)) ck++;
      if (naArr[j] >= oa) ca++;
    }
    const pKS_emp = (ck + 1) / (nPerms + 1);
    const pAD_emp = (ca + 1) / (nPerms + 1);

    // NES
    const nullKSarr = Array.from(nkArr);
    const nullADarr = Array.from(naArr);
    const nes    = calcNES(ok, nullKSarr);
    const nes_ad = calcNES_AD(oa, nullADarr);

    let pKS_fit = pKS_emp, pAD_fit = pAD_emp;
    let gammaFit = null, ggFit = null;

    if (engine === 'gg') {
      // Fit Gamma to |null KS|
      const absNullKS = nkArr.map(v => Math.abs(v));
      gammaFit = fitGamma(Array.from(absNullKS));
      pKS_fit  = pvalFromGamma(ok, gammaFit, pKS_emp);

      // Fit Generalised Gamma to null AD
      ggFit   = fitGenGamma(nullADarr);
      pAD_fit = pvalFromGG(oa, ggFit, pAD_emp);
    }

    const pCauchy = cauchyCombine(pKS_fit, pAD_fit);

    results.push({
      name:     p.name,
      url:      p.url || null,
      size:     p.size,
      es:       ok,
      ad:       oa,
      nes,
      nes_ad,
      pKS_emp, pAD_emp,
      pKS_fit, pAD_fit,
      pCauchy,
      gammaFit, ggFit,
      curve:    obsStats[pi].curve,
      obsOrd:   Array.from(obsOrd),
      peakIdx:  obsStats[pi].peakIdx
    });
  }

  // ── 4. FDR (BH) on primary p-value (pCauchy) ──────────
  if (results.length > 1) {
    const pvec  = results.map(r => r.pCauchy);
    const qvec  = bhFDR(pvec);
    results.forEach((r, i) => { r.fdr = qvec[i]; });
  }

  onProgress(100, `Done — ${results.length} pathway(s) analysed`);
  return results;
}

function yieldToBrowser() {
  return new Promise(r => setTimeout(r, 0));
}
