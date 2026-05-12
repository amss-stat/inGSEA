// ═══════════════════════════════════════════════════════════
//  core/gsea.js  ·  v2.4
//
//  Changes:
//  • Passes weight exponent p to calcGSEAStats
//  • permSNR as Float32Array for faster sort (ranking only)
//  • ord as plain Array for V8 TimSort fast path
//  • Time-adaptive chunk sizing
//  • Abort token
// ═══════════════════════════════════════════════════════════
'use strict';

import {
  calcSNR, rankOrder, calcGSEAStats,
  cauchyCombine, bhFDR, shuffle,
  calcNES, calcNES_AD, empP_KS, empP_AD
} from './math.js';

import { pvalGG, pvalGamma } from './distributions.js';

const MAX_PERMS       = 2000;
const TARGET_CHUNK_MS = 40;

export async function runGSEA(opts) {
  const { exprMat, geneNames, pathways, engine, onProgress, abortSignal } = opts;
  const nPerms  = Math.min(opts.nPerms | 0, MAX_PERMS);
  const nG      = geneNames.length;
  const nS      = exprMat[0].length;
  const nCase   = opts.nCase | 0;
  const nP      = pathways.length;
  const wt      = typeof opts.weightP === 'number' ? opts.weightP : 1;

  if (nCase < 2 || nCase > nS - 2)
    throw new Error(`nCase=${nCase} out of range [2, ${nS-2}] (nSamples=${nS})`);

  const _ab = () => abortSignal?.aborted === true;

  const cIdx = new Uint16Array(nCase);
  const tIdx = new Uint16Array(nS - nCase);
  for (let i = 0; i < nCase; i++)      cIdx[i] = i;
  for (let i = 0; i < nS - nCase; i++) tIdx[i] = i + nCase;

  // ── 1. Observed ──────────────────────────────────────────────
  onProgress(0, 'Observed', 'Computing SNR & enrichment…');
  // Observed SNR uses Float64 for full precision (stored for NES)
  const obsSNR = new Float64Array(nG);
  calcSNR(exprMat, nG, cIdx, tIdx, obsSNR);

  // ord as plain Array — V8 TimSort is faster than TypedArray.sort
  const obsOrd = new Array(nG);
  rankOrder(obsSNR, nG, obsOrd);

  const obsStats = pathways.map(p =>
    calcGSEAStats(obsSNR, obsOrd, p.mask, nG, true, wt)
  );

  if (_ab()) throw new Error('Aborted');

  // ── 2. Permutations ──────────────────────────────────────────
  const nullKS = Array.from({ length: nP }, () => new Float64Array(nPerms));
  const nullAD = Array.from({ length: nP }, () => new Float64Array(nPerms));

  const permArr = new Uint16Array(nS);
  for (let i = 0; i < nS; i++) permArr[i] = i;

  // Float32Array for permutation SNR: halves memory bandwidth
  // for the sort step (ranking only, precision sufficient)
  const permSNR = new Float32Array(nG);

  // Plain Array for ord (V8 TimSort fast path)
  const permOrd = new Array(nG);

  const permC = permArr.subarray(0, nCase);   // zero-copy views
  const permT = permArr.subarray(nCase);

  const t0 = performance.now();
  let chunkSize = 5;
  let done = 0;

  while (done < nPerms) {
    if (_ab()) throw new Error('Aborted');

    const end      = Math.min(done + chunkSize, nPerms);
    const tChunk   = performance.now();

    for (let j = done; j < end; j++) {
      shuffle(permArr);
      // Write into permSNR in-place (Float32, no alloc)
      calcSNR(exprMat, nG, permC, permT, permSNR);
      rankOrder(permSNR, nG, permOrd);

      for (let pi = 0; pi < nP; pi++) {
        const st = calcGSEAStats(permSNR, permOrd, pathways[pi].mask, nG, false, wt);
        nullKS[pi][j] = st.ks;
        nullAD[pi][j] = st.ad;
      }
    }

    done = end;

    // Adapt chunk size to TARGET_CHUNK_MS
    const elapsed = performance.now() - tChunk;
    if (elapsed > 0 && done < nPerms) {
      const newChunk = Math.round(chunkSize * TARGET_CHUNK_MS / elapsed);
      chunkSize = Math.max(1, Math.min(500, newChunk));
    }

    const totalElapsed = (performance.now() - t0) / 1000;
    const pct   = Math.round(done / nPerms * 80);
    const eta   = totalElapsed > 0.1
      ? ((nPerms - done) / (done / totalElapsed)).toFixed(0) : '—';
    const rate  = (done / (totalElapsed + 0.001)).toFixed(0);
    onProgress(pct, 'Permutations',
      `${done}/${nPerms} · ETA ${eta}s · ${rate} perm/s`);

    await _yield();
  }

  if (_ab()) throw new Error('Aborted');

  // ── 3. Empirical p & NES ─────────────────────────────────────
  onProgress(81, 'Statistics', 'p-values & NES…');
  await _yield();

  const empArr = pathways.map((p, pi) => {
    const ok = obsStats[pi].ks;
    const oa = obsStats[pi].ad;
    return {
      ok, oa,
      pKS_emp: empP_KS(ok, nullKS[pi], nPerms),
      pAD_emp: empP_AD(oa, nullAD[pi], nPerms),
      nes:     calcNES(ok, nullKS[pi]),
      nes_ad:  calcNES_AD(oa, nullAD[pi])
    };
  });

  if (_ab()) throw new Error('Aborted');

  // ── 4. Parametric fitting ─────────────────────────────────────
  let pKS_arr, pAD_arr;
  const useParam = engine === 'parametric';

  if (useParam) {
    onProgress(85, 'Fitting', 'Fitting null distributions (jStat)…');
    await _yield();
    pKS_arr = [];
    pAD_arr = [];
    for (let pi = 0; pi < nP; pi++) {
      if (_ab()) throw new Error('Aborted');
      const r = empArr[pi];
      pKS_arr.push(pvalGamma(r.ok, nullKS[pi], r.pKS_emp));
      pAD_arr.push(pvalGG(r.oa, nullAD[pi], r.pAD_emp));
      if (pi % 3 === 2) await _yield();
    }
  } else {
    pKS_arr = empArr.map(r => r.pKS_emp);
    pAD_arr = empArr.map(r => r.pAD_emp);
  }

  // ── 5. Assemble ──────────────────────────────────────────────
  onProgress(97, 'Assembling', 'Cauchy & FDR…');
  await _yield();

  const results = pathways.map((p, pi) => {
    const e  = empArr[pi];
    const pC = cauchyCombine(pKS_arr[pi], pAD_arr[pi]);
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
      pKS:     pKS_arr[pi],
      pAD:     pAD_arr[pi],
      pCauchy: pC,
      fdr:     null,
      curve:   obsStats[pi].curve,
      obsOrd,
      peakIdx: obsStats[pi].peakIdx
    };
  });

  // BH-FDR on pCauchy
  if (results.length > 1) {
    const pv = results.map(r => r.pCauchy);
    const qv = bhFDR(pv);
    results.forEach((r, i) => { r.fdr = qv[i]; });
  }

  onProgress(100, 'Done', `${results.length} pathway(s) complete`);
  return results;
}

const _yield = () => new Promise(r => setTimeout(r, 0));
