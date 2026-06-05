'use strict';

import {
  calcSNR, rankOrder, calcGSEAStats,
  cauchyCombine, shuffle,
  calcNES, calcNES_AD, empP_KS, empP_AD,
  calcGseaFDR
} from './math.js';

import { pvalGG } from './distributions.js';

const MAX_PERMS       = 5000;
const TARGET_CHUNK_MS = 16;

export async function runGSEA(opts) {
  const {
    exprMat, geneNames, pathways,
    engine, onProgress, abortSignal
  } = opts;

  const nPerms = Math.min(opts.nPerms | 0, MAX_PERMS);
  const nG     = geneNames.length;
  const nS     = exprMat[0].length;
  const nCase  = opts.nCase | 0;
  const nP     = pathways.length;
  const wt     = typeof opts.weightP === 'number' ? opts.weightP : 1;

  if (nCase < 2 || nCase > nS - 2)
    throw new Error(
      `nCase=${nCase} out of range [2, ${nS - 2}] (nSamples=${nS})`
    );

  const _ab = () => abortSignal?.aborted === true;

  const cIdx = new Uint16Array(nCase);
  const tIdx = new Uint16Array(nS - nCase);
  for (let i = 0; i < nCase; i++)      cIdx[i] = i;
  for (let i = 0; i < nS - nCase; i++) tIdx[i] = i + nCase;

  // ── 1. Observed ───────────────────────────────────────────
  onProgress(0, 'Observed', 'Computing SNR & enrichment…');

  const obsSNR = new Float64Array(nG);
  calcSNR(exprMat, nG, cIdx, tIdx, obsSNR);

  const obsOrd = new Array(nG);
  rankOrder(obsSNR, nG, obsOrd);

  const obsStats = pathways.map(p =>
    calcGSEAStats(obsSNR, obsOrd, p.mask, nG, true, wt)
  );

  if (_ab()) throw new Error('Aborted');

  // ── 2. Permutations ───────────────────────────────────────
  const nullKS = Array.from({ length: nP }, () => new Float64Array(nPerms));
  const nullAD = Array.from({ length: nP }, () => new Float64Array(nPerms));

  const permArr = new Uint16Array(nS);
  for (let i = 0; i < nS; i++) permArr[i] = i;

  const permSNR = new Float32Array(nG);
  const permOrd = new Array(nG);
  const permC   = permArr.subarray(0, nCase);
  const permT   = permArr.subarray(nCase);

  const t0 = performance.now();
  let chunkSize = 5;
  let done = 0;

  while (done < nPerms) {
    if (_ab()) throw new Error('Aborted');

    const end    = Math.min(done + chunkSize, nPerms);
    const tChunk = performance.now();

    for (let j = done; j < end; j++) {
      shuffle(permArr);
      calcSNR(exprMat, nG, permC, permT, permSNR);
      rankOrder(permSNR, nG, permOrd);
      for (let pi = 0; pi < nP; pi++) {
        const st = calcGSEAStats(
          permSNR, permOrd, pathways[pi].mask, nG, false, wt
        );
        nullKS[pi][j] = st.ks;
        nullAD[pi][j] = st.ad;
      }
    }

    done = end;

    const elapsed = performance.now() - tChunk;
    if (elapsed > 0 && done < nPerms) {
      const newChunk = Math.round(chunkSize * TARGET_CHUNK_MS / elapsed);
      chunkSize = Math.max(1, Math.min(500, newChunk));
    }

    const totalElapsed = (performance.now() - t0) / 1000;
    const pct  = Math.round(done / nPerms * 75);
    onProgress(pct, 'Permutations', `${done}/${nPerms} permutations`);

    await _yield();
  }

  if (_ab()) throw new Error('Aborted');

  // ── 3. Empirical p-values & NES ───────────────────────────
  onProgress(76, 'Statistics', 'Empirical p-values & NES…');
  await _yield();

  const empArr = new Array(nP);
  for (let pi = 0; pi < nP; pi++) {
    if (_ab()) throw new Error('Aborted');

    const ok = obsStats[pi].ks;
    const oa = obsStats[pi].ad;
    empArr[pi] = {
      ok,
      oa,
      pKS_emp: empP_KS(ok, nullKS[pi], nPerms),
      pAD_emp: empP_AD(oa, nullAD[pi], nPerms),
      nes:     calcNES(ok, nullKS[pi]),
      nes_ad:  calcNES_AD(oa, nullAD[pi])
    };
  }

  if (_ab()) throw new Error('Aborted');

  // ── 4. AD parametric fitting (engine = 'parametric') ──────
  const useParam = engine === 'parametric';
  const pKS_arr  = new Array(nP);
  const pAD_arr  = new Array(nP);

  for (let pi = 0; pi < nP; pi++) {
    pKS_arr[pi] = empArr[pi].pKS_emp;
  }

  if (useParam) {
    onProgress(80, 'Fitting', 'Fitting AD null distributions (GG)…');
    await _yield();
    for (let pi = 0; pi < nP; pi++) {
      if (_ab()) throw new Error('Aborted');
      const r   = empArr[pi];
      const res = pvalGG(r.oa, nullAD[pi], r.pAD_emp);
      pAD_arr[pi] = res.fitted ? res.p : r.pAD_emp;
      if (pi % 10 === 0) await _yield();
    }
  } else {
    for (let pi = 0; pi < nP; pi++) pAD_arr[pi] = empArr[pi].pAD_emp;
  }

  // ── 5. Cauchy combination ─────────────────────────────────
  const pCauchy_arr = new Array(nP);
  for (let pi = 0; pi < nP; pi++) {
    pCauchy_arr[pi] = cauchyCombine(pKS_arr[pi], pAD_arr[pi]);
  }

  // ── 6. GSEA FDR ───────────────────────────────────────────
  onProgress(93, 'FDR', 'Building Global NES Distributions…');
  await _yield();

  const obsNES_ks = new Float64Array(nP);
  const obsNES_ad = new Float64Array(nP);
  const nullNES_ks = Array.from({ length: nP }, () => new Float64Array(nPerms));
  const nullNES_ad = Array.from({ length: nP }, () => new Float64Array(nPerms));

  for (let pi = 0; pi < nP; pi++) {
    const e = empArr[pi];
    obsNES_ks[pi] = e.nes;
    obsNES_ad[pi] = e.nes_ad;

    let sumPosKS = 0, countPosKS = 0;
    let sumNegKS = 0, countNegKS = 0;
    let sumAD = 0;

    for (let j = 0; j < nPerms; j++) {
      const ks = nullKS[pi][j];
      if (ks >= 0) { sumPosKS += ks; countPosKS++; }
      else         { sumNegKS += ks; countNegKS++; }
      sumAD += nullAD[pi][j];
    }

    const meanPosKS = countPosKS > 0 ? sumPosKS / countPosKS : 0;
    const meanNegKS = countNegKS > 0 ? Math.abs(sumNegKS / countNegKS) : 0;
    const meanAD    = sumAD / nPerms;

    for (let j = 0; j < nPerms; j++) {
      const ks = nullKS[pi][j];
      nullNES_ks[pi][j] = ks >= 0 
        ? (meanPosKS > 0 ? ks / meanPosKS : 0)
        : (meanNegKS > 0 ? ks / meanNegKS : 0);
      nullNES_ad[pi][j] = meanAD > 0 ? nullAD[pi][j] / meanAD : 0;
    }
  }

  onProgress(97, 'FDR', 'Computing Empirical FDRs…');
  await _yield();

  const fdr_ks_arr = calcGseaFDR(obsNES_ks, nullNES_ks, true);
  const fdr_ad_arr = calcGseaFDR(obsNES_ad, nullNES_ad, false);

  const results = pathways.map((p, pi) => {
    const e = empArr[pi];
    return {
      name:    p.name,
      url:     p.url ?? null,
      size:    p.size,

      es:      e.ok,
      ad:      e.oa,
      nes:     e.nes,
      nes_ad:  e.nes_ad,

      pKS:     pKS_arr[pi],
      pAD:     pAD_arr[pi],
      pCauchy: pCauchy_arr[pi],

      fdr_ks:  fdr_ks_arr[pi],
      fdr_ad:  fdr_ad_arr[pi],

      curve:   obsStats[pi].curve,
      obsOrd,
      peakIdx: obsStats[pi].peakIdx
    };
  });

  onProgress(100, 'Done', `${results.length} pathway(s) complete`);
  return results;
}

const _yield = () => new Promise(r => setTimeout(r, 0));
