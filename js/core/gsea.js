// ═══════════════════════════════════════════════════════════
//  core/gsea.js  ·  v2.6
//  Key fix: pKS_par / pAD_par stored separately from
//  pKS_emp / pAD_emp so CSV and table can show both.
//  pKS / pAD = parametric when engine=parametric, else empirical.
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

  // Convert Float64Array null distributions to plain Arrays
  // once here — avoids repeated conversion in fitting loops
  const nullKS_arr = nullKS.map(a => Array.from(a));
  const nullAD_arr = nullAD.map(a => Array.from(a));

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

  // ── 4. Parametric fitting ─────────────────────────────────
  const useParam = engine === 'parametric';
  // pKS_arr / pAD_arr hold the p-values used for Cauchy combination:
  //   parametric when engine='parametric' AND fit succeeds
  //   empirical  as fallback (or always when engine='permutation')
  const pKS_arr = new Array(nP);
  const pAD_arr = new Array(nP);
  // pKS_par_arr / pAD_par_arr: null when engine='permutation'
  const pKS_par_arr = new Array(nP).fill(null);
  const pAD_par_arr = new Array(nP).fill(null);

  if (useParam) {
    onProgress(85, 'Fitting', 'Fitting null distributions…');
    await _yield();

    for (let pi = 0; pi < nP; pi++) {
      if (_ab()) throw new Error('Aborted');

      const r = empArr[pi];

      // pvalGamma / pvalGG return the parametric p-value,
      // or fall back to empP if fitting fails.
      // We capture both outcomes separately.
      const pKS_p = pvalGamma(r.ok, nullKS_arr[pi], r.pKS_emp);
      const pAD_p = pvalGG(r.oa,   nullAD_arr[pi], r.pAD_emp);

      // Detect fallback: if returned value === empP exactly,
      // the parametric fit failed and we fell back.
      // Store null in _par fields to make this explicit in output.
      pKS_par_arr[pi] = (pKS_p === r.pKS_emp) ? null : pKS_p;
      pAD_par_arr[pi] = (pAD_p === r.pAD_emp) ? null : pAD_p;

      // For Cauchy: use parametric if available, else empirical
      pKS_arr[pi] = pKS_p;
      pAD_arr[pi] = pAD_p;

      if (pi % 3 === 2) await _yield();
    }
  } else {
    for (let pi = 0; pi < nP; pi++) {
      pKS_arr[pi] = empArr[pi].pKS_emp;
      pAD_arr[pi] = empArr[pi].pAD_emp;
      // pKS_par_arr / pAD_par_arr remain null
    }
  }

  // ── 5. Assemble results ───────────────────────────────────
  onProgress(97, 'Assembling', 'Cauchy & FDR…');
  await _yield();

  const results = pathways.map((p, pi) => {
    const e  = empArr[pi];
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
      // Parametric — null when engine=permutation or fit failed
      pKS_par:  pKS_par_arr[pi],
      pAD_par:  pAD_par_arr[pi],
      // Combined — used for Cauchy and displayed as primary p-value
      pKS:      pKS_arr[pi],
      pAD:      pAD_arr[pi],
      pCauchy:  pC,
      fdr:      null,   // filled below
      curve:    obsStats[pi].curve,
      obsOrd,
      peakIdx:  obsStats[pi].peakIdx
    };
  });

  // ── 6. BH-FDR on pCauchy ─────────────────────────────────
  // Always compute FDR when ≥ 2 pathways; UI decides whether to show it
  if (results.length >= 2) {
    const pv = results.map(r => r.pCauchy);
    const qv = bhFDR(pv);
    results.forEach((r, i) => { r.fdr = qv[i]; });
  }

  onProgress(100, 'Done', `${results.length} pathway(s) complete`);
  return results;
}

const _yield = () => new Promise(r => setTimeout(r, 0));
