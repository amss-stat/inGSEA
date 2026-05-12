// ═══════════════════════════════════════════════════════════
//  core/gsea.js  ·  v2.7
//
//  Key changes from v2.6:
//  • KS always uses empirical p-values (no parametric fitting)
//  • AD uses GG parametric approximation when engine='parametric'
//  • pvalGG returns { p, fitted } — no sentinel-value detection
//  • Cauchy combines pKS_emp with pAD (parametric or empirical)
//  • Removed pvalGamma import (no longer used)
//  • pKS_par always null (KS has no parametric path)
// ═══════════════════════════════════════════════════════════
'use strict';

import {
  calcSNR, rankOrder, calcGSEAStats,
  cauchyCombine, bhFDR, shuffle,
  calcNES, calcNES_AD, empP_KS, empP_AD
} from './math.js';

import { pvalGG } from './distributions.js';

const MAX_PERMS       = 2000;
const TARGET_CHUNK_MS = 40;

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
    const pct  = Math.round(done / nPerms * 80);
    const eta  = totalElapsed > 0.1
      ? ((nPerms - done) / (done / totalElapsed)).toFixed(0) : '—';
    const rate = (done / (totalElapsed + 0.001)).toFixed(0);
    onProgress(pct, 'Permutations',
      `${done}/${nPerms} · ETA ${eta}s · ${rate} perm/s`);

    await _yield();
  }

  if (_ab()) throw new Error('Aborted');

  // ── 3. Empirical p & NES ──────────────────────────────────
  onProgress(81, 'Statistics', 'p-values & NES…');
  await _yield();

  const empArr = pathways.map((p, pi) => {
    const ok = obsStats[pi].ks;
    const oa = obsStats[pi].ad;
    return {
      ok,
      oa,
      pKS_emp: empP_KS(ok, nullKS[pi], nPerms),
      pAD_emp: empP_AD(oa, nullAD[pi], nPerms),
      nes:     calcNES(ok, nullKS[pi]),
      nes_ad:  calcNES_AD(oa, nullAD[pi])
    };
  });

  if (_ab()) throw new Error('Aborted');

  // ── 4. AD parametric fitting (when engine = 'parametric') ─
  const useParam = engine === 'parametric';

  // KS: always empirical — no parametric path
  const pKS_arr = new Array(nP);
  for (let pi = 0; pi < nP; pi++) {
    pKS_arr[pi] = empArr[pi].pKS_emp;
  }

  // AD: parametric GG when requested, empirical as fallback
  const pAD_arr     = new Array(nP);
  const pAD_par_arr = new Array(nP).fill(null);

  if (useParam) {
    onProgress(85, 'Fitting', 'Fitting AD null distributions (GG)…');
    await _yield();

    // Convert Float64Array null distributions to plain Arrays once
    const nullAD_arr = nullAD.map(a => Array.from(a));

    let fitSuccess = 0;
    for (let pi = 0; pi < nP; pi++) {
      if (_ab()) throw new Error('Aborted');

      const r   = empArr[pi];
      const res = pvalGG(r.oa, nullAD_arr[pi], r.pAD_emp);

      if (res.fitted) {
        pAD_par_arr[pi] = res.p;
        pAD_arr[pi]     = res.p;
        fitSuccess++;
      } else {
        pAD_par_arr[pi] = null;
        pAD_arr[pi]     = r.pAD_emp;
      }

      if (pi % 3 === 2) await _yield();
    }

    console.log(`GG fit success: ${fitSuccess}/${nP} pathways`);
  } else {
    for (let pi = 0; pi < nP; pi++) {
      pAD_arr[pi] = empArr[pi].pAD_emp;
    }
  }

  // ── 5. Assemble results ───────────────────────────────────
  onProgress(97, 'Assembling', 'Cauchy & FDR…');
  await _yield();

  const results = pathways.map((p, pi) => {
    const e  = empArr[pi];
    // Cauchy combines: pKS (always empirical) + pAD (parametric or empirical)
    const pC = cauchyCombine(pKS_arr[pi], pAD_arr[pi]);
    return {
      name:     p.name,
      url:      p.url ?? null,
      size:     p.size,
      es:       e.ok,
      ad:       e.oa,
      nes:      e.nes,
      nes_ad:   e.nes_ad,
      // Empirical — always present
      pKS_emp:  e.pKS_emp,
      pAD_emp:  e.pAD_emp,
      // Parametric — KS is always null; AD is null when not fitted
      pKS_par:  null,
      pAD_par:  pAD_par_arr[pi],
      // Primary p-values used for display and Cauchy
      pKS:      pKS_arr[pi],      // always empirical
      pAD:      pAD_arr[pi],      // parametric if fitted, else empirical
      pCauchy:  pC,
      fdr:      null,              // filled below
      curve:    obsStats[pi].curve,
      obsOrd,
      peakIdx:  obsStats[pi].peakIdx
    };
  });

  // ── 6. BH-FDR on pCauchy ─────────────────────────────────
  if (results.length >= 2) {
    const pv = results.map(r => r.pCauchy);
    const qv = bhFDR(pv);
    results.forEach((r, i) => { r.fdr = qv[i]; });
  }

  onProgress(100, 'Done', `${results.length} pathway(s) complete`);
  return results;
}

const _yield = () => new Promise(r => setTimeout(r, 0));
