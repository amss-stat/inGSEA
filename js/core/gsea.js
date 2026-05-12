// ═══════════════════════════════════════════════════════════
//  core/gsea.js  ·  v2.9
//
//  Pipeline:
//    1. Observed SNR + enrichment stats
//    2. Permutations → nullKS[nP][nPerms], nullAD[nP][nPerms]
//    3. Empirical p-values (KS two-sided, AD one-sided)
//    4. AD parametric fitting via GG (engine='parametric')
//    5. Cauchy combination → pCauchy
//    6. Assemble results
//    7. FDR: Storey q-value on pCauchy
//
//  KS p-value: always empirical (KS null ES is signed,
//              gamma approximation inapplicable)
//  AD p-value: parametric GG when engine='parametric',
//              empirical fallback when fit fails
//  FDR:        Storey q-value on pCauchy; estimates π₀
//              (proportion of true nulls) for improved power
//              over BH under positive correlation
// ═══════════════════════════════════════════════════════════
'use strict';

import {
  calcSNR, rankOrder, calcGSEAStats,
  cauchyCombine, shuffle,
  calcNES, calcNES_AD, empP_KS, empP_AD,
  storeyQ
} from './math.js';

import { pvalGG } from './distributions.js';

const MAX_PERMS       = 2000;
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
  // nullKS[pi][j] = KS enrichment score for pathway pi,
  //                 permutation j  (signed, can be negative)
  // nullAD[pi][j] = AD statistic  (always non-negative)
  const nullKS = Array.from({ length: nP }, () => new Float64Array(nPerms));
  const nullAD = Array.from({ length: nP }, () => new Float64Array(nPerms));

  const permArr = new Uint16Array(nS);
  for (let i = 0; i < nS; i++) permArr[i] = i;

  // Float32 for permutation SNR: halves memory bandwidth,
  // precision sufficient for ranking only
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
    const eta  = totalElapsed > 0.1
      ? ((nPerms - done) / (done / totalElapsed)).toFixed(0) : '—';
    const rate = (done / (totalElapsed + 0.001)).toFixed(0);
    onProgress(pct, 'Permutations',
      `${done}/${nPerms} · ETA ${eta}s · ${rate} perm/s`);

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

    if (pi % 50 === 49) {
      onProgress(76, 'Statistics',
        `Empirical p-values (${pi + 1}/${nP})…`);
      await _yield();
    }
  }

  if (_ab()) throw new Error('Aborted');

  // ── 4. AD parametric fitting (engine = 'parametric') ──────
  // KS always uses empirical p-values:
  //   null KS ES is signed (positive and negative), making
  //   gamma approximation inapplicable.
  // AD uses GG parametric approximation when requested:
  //   AD statistic is always non-negative, GG fit is appropriate.
  //   Falls back to empirical if fit fails.
  const useParam    = engine === 'parametric';
  const pKS_arr     = new Array(nP);
  const pAD_arr     = new Array(nP);
  const pAD_par_arr = new Array(nP).fill(null);

  // KS: always empirical
  for (let pi = 0; pi < nP; pi++) {
    pKS_arr[pi] = empArr[pi].pKS_emp;
  }

  if (useParam) {
    onProgress(80, 'Fitting', 'Fitting AD null distributions (GG)…');
    await _yield();

    let fitSuccess = 0;
    for (let pi = 0; pi < nP; pi++) {
      if (_ab()) throw new Error('Aborted');

      const r   = empArr[pi];
      const res = pvalGG(r.oa, nullAD[pi], r.pAD_emp);

      if (res.fitted) {
        pAD_par_arr[pi] = res.p;
        pAD_arr[pi]     = res.p;
        fitSuccess++;
      } else {
        pAD_arr[pi] = r.pAD_emp;
      }

      // Yield every 5 pathways: each pvalGG call runs Nelder-Mead
      // (up to 1000 iterations), so 5 pathways ≈ 5-25ms
      if (pi % 5 === 4 || pi === nP - 1) {
        const pct = 80 + Math.round(((pi + 1) / nP) * 10);
        onProgress(pct, 'Fitting',
          `AD fitting (${pi + 1}/${nP}) · fitted: ${fitSuccess}`);
        await _yield();
      }
    }

    console.log(`GG fit: ${fitSuccess}/${nP} pathways fitted parametrically`);
  } else {
    for (let pi = 0; pi < nP; pi++) {
      pAD_arr[pi] = empArr[pi].pAD_emp;
    }
  }

  if (_ab()) throw new Error('Aborted');

  // ── 5. Cauchy combination ─────────────────────────────────
  // Combines pKS (empirical) and pAD (parametric or empirical)
  // into a single p-value per pathway.
  // Liu & Xie (2020): robust to correlation between input p-values.
  onProgress(91, 'Combining', 'Cauchy combination…');
  await _yield();

  const pCauchy_arr = new Array(nP);
  for (let pi = 0; pi < nP; pi++) {
    pCauchy_arr[pi] = cauchyCombine(pKS_arr[pi], pAD_arr[pi]);
  }

  // ── 6. Assemble results ───────────────────────────────────
  // FDR is filled in step 7 after assembly.
  onProgress(93, 'Assembling', 'Assembling results…');
  await _yield();

  const results = pathways.map((p, pi) => {
    const e = empArr[pi];
    return {
      name:    p.name,
      url:     p.url ?? null,
      size:    p.size,

      // Effect sizes
      es:      e.ok,       // raw enrichment score (signed)
      ad:      e.oa,       // Anderson-Darling statistic
      nes:     e.nes,      // normalised ES (KS)
      nes_ad:  e.nes_ad,   // normalised AD

      // Primary p-values (displayed in main table columns)
      // pKS: always empirical permutation p-value
      // pAD: parametric GG if fitted, else empirical
      pKS:     pKS_arr[pi],
      pAD:     pAD_arr[pi],
      pCauchy: pCauchy_arr[pi],

      // Empirical p-values (always available, extended columns)
      pKS_emp: e.pKS_emp,
      pAD_emp: e.pAD_emp,

      // Parametric AD p-value: null when engine='permutation'
      // or when GG fit failed (extended columns)
      pKS_par: null,              // KS has no parametric path
      pAD_par: pAD_par_arr[pi],

      // FDR: filled below after assembly
      fdr:     null,

      // Curve data (observed enrichment walk)
      curve:   obsStats[pi].curve,
      obsOrd,
      peakIdx: obsStats[pi].peakIdx
    };
  });

  // ── 7. FDR (Storey q-value on pCauchy) ───────────────────
  // storeyQ estimates π₀ (proportion of true nulls) from the
  // pCauchy distribution, then applies π₀-weighted BH.
  // More powerful than plain BH under positive correlation;
  // reduces to BH when π₀ cannot be estimated reliably (nP < 10).
  onProgress(97, 'FDR', 'Storey q-value…');
  await _yield();

  if (nP >= 2) {
    const qArr = storeyQ(pCauchy_arr);
    results.forEach((r, i) => { r.fdr = qArr[i]; });
  } else {
    // Single pathway: FDR not meaningful, mirror pCauchy
    results[0].fdr = results[0].pCauchy;
  }

  onProgress(100, 'Done', `${results.length} pathway(s) complete`);
  return results;
}

const _yield = () => new Promise(r => setTimeout(r, 0));
