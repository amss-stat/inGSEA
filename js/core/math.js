// ═══════════════════════════════════════════════════════════
//  core/math.js  ·  v2.4
//  Changes:
//  • calcSNR: accepts optional Float32Array output (faster sort)
//  • calcGSEAStats: accepts weight exponent p (default 1)
//  • bhFDR: clamp q to [0,1]; correct sweep
// ═══════════════════════════════════════════════════════════
'use strict';

/**
 * Welch SNR: (μ_case − μ_ctrl) / (σ_case + σ_ctrl + ε)
 *
 * Performance notes:
 * - out can be Float64Array or Float32Array.
 *   Float32 halves memory bandwidth for the sort step.
 *   Precision is sufficient for ranking (we only need the order).
 * - cIdx / tIdx may be Uint16Array views (subarray) — works correctly.
 */
export function calcSNR(mat, nG, cIdx, tIdx, out) {
  const snr = out ?? new Float64Array(nG);
  const nc = cIdx.length, nt = tIdx.length;
  for (let g = 0; g < nG; g++) {
    const row = mat[g];
    let m1 = 0, m2 = 0;
    for (let i = 0; i < nc; i++) m1 += row[cIdx[i]];
    for (let i = 0; i < nt; i++) m2 += row[tIdx[i]];
    m1 /= nc; m2 /= nt;
    let v1 = 0, v2 = 0;
    for (let i = 0; i < nc; i++) { const d = row[cIdx[i]] - m1; v1 += d * d; }
    for (let i = 0; i < nt; i++) { const d = row[tIdx[i]] - m2; v2 += d * d; }
    const s1 = nc > 1 ? Math.sqrt(v1 / (nc - 1)) : 0;
    const s2 = nt > 1 ? Math.sqrt(v2 / (nt - 1)) : 0;
    snr[g] = (m1 - m2) / (s1 + s2 + 1e-9);
  }
  return snr;
}

/**
 * Descending argsort of snr into pre-allocated Array ord.
 * Uses plain Array (not TypedArray) for ord so V8 uses its
 * optimised TimSort path rather than the TypedArray fallback.
 *
 * @param {Float64Array|Float32Array} snr
 * @param {number}  nG
 * @param {number[]} ord  plain Array, length nG, pre-allocated
 */
export function rankOrder(snr, nG, ord) {
  for (let i = 0; i < nG; i++) ord[i] = i;
  // V8 TimSort on plain Array with numeric comparator is faster
  // than Int32Array.sort with a closure in most engines.
  ord.sort((a, b) => snr[b] - snr[a]);
}

/**
 * GSEA enrichment walk — standard weighted KS statistic.
 *
 * Weight exponent p (default 1):
 *   p = 0 : unweighted (classical KS, equal weight to all hits)
 *   p = 1 : weighted by |SNR|  (original GSEA default)
 *   p = 2 : quadratic weighting
 *
 * For each position i in the ranked list:
 *   hit  gene: Δ = |snr_i|^p / N_R   where N_R = Σ_{j∈S} |snr_j|^p
 *   miss gene: Δ = -1 / (N - N_H)    where N_H = |S|, N = total genes
 *
 * ES = running sum value at its maximum absolute deviation.
 * AD = Σ_{i=0}^{N-2}  ES_i² / (F_i · (1-F_i))
 *   where F_i = (i+1)/N  (fraction of genes processed through step i+1)
 *
 * @param {Float64Array|Float32Array} snr    ranked-gene SNR values
 * @param {number[]|Int32Array}       ord    indices sorted by descending SNR
 * @param {Uint8Array}                mask   1 = gene in set
 * @param {number}                    nG
 * @param {boolean}                   wantCurve
 * @param {number}                    [p=1]  weight exponent
 * @returns {{ ks, ad, curve, peakIdx }}
 */
export function calcGSEAStats(snr, ord, mask, nG, wantCurve, p = 1) {
  // Pre-compute hit count and weighted sum N_R
  let nHits = 0, hitSum = 0;
  for (let i = 0; i < nG; i++) {
    if (mask[ord[i]]) {
      nHits++;
      const absV = snr[ord[i]] < 0 ? -snr[ord[i]] : snr[ord[i]];
      hitSum += p === 1 ? absV
              : p === 0 ? 1
              : Math.pow(absV, p);
    }
  }

  if (nHits === 0 || nHits === nG) {
    return { ks: 0, ad: 0, peakIdx: 0,
             curve: wantCurve ? new Float64Array(nG) : null };
  }

  const nMiss    = nG - nHits;
  const missStep = -1.0 / nMiss;
  const invHit   = 1.0 / (hitSum + 1e-12);
  const curve    = wantCurve ? new Float64Array(nG) : null;

  let cum = 0.0, ks = 0.0, maxAbs = 0.0, peakIdx = 0, ad = 0.0;

  for (let i = 0; i < nG; i++) {
    if (mask[ord[i]]) {
      const absV = snr[ord[i]] < 0 ? -snr[ord[i]] : snr[ord[i]];
      cum += (p === 1 ? absV : p === 0 ? 1 : Math.pow(absV, p)) * invHit;
    } else {
      cum += missStep;
    }

    if (wantCurve) curve[i] = cum;

    const a = cum < 0 ? -cum : cum;
    if (a > maxAbs) { maxAbs = a; ks = cum; peakIdx = i; }

    // AD: interior positions only (i = 0 … nG-2)
    if (i < nG - 1) {
      const F = (i + 1) / nG;
      ad += (cum * cum) / (F * (1.0 - F));
    }
  }

  return { ks, ad, curve, peakIdx };
}

/**
 * Cauchy combination of p-values (equal weights).
 * Liu & Xie 2020.
 */
export function cauchyCombine(...pvals) {
  let s = 0;
  for (let p of pvals) {
    p = p < 1e-16 ? 1e-16 : p > 1 - 1e-16 ? 1 - 1e-16 : p;
    s += Math.tan((0.5 - p) * Math.PI);
  }
  s /= pvals.length;
  return 0.5 - Math.atan(s) / Math.PI;
}

/**
 * Benjamini–Hochberg FDR (step-up procedure).
 * Returns Float64Array of q-values aligned to input order.
 * All q-values are clamped to [0, 1].
 *
 * Algorithm:
 *   Sort p ascending → p(1) ≤ p(2) ≤ … ≤ p(m)
 *   BH q(i) = min_{j≥i} { p(j) · m / j }   (cumulative min from right)
 */
export function bhFDR(pvals) {
  const n   = pvals.length;
  if (n === 0) return new Float64Array(0);
  if (n === 1) return Float64Array.from([Math.min(pvals[0], 1)]);

  const idx = Array.from({ length: n }, (_, i) => i);
  idx.sort((a, b) => pvals[a] - pvals[b]);

  const q = new Float64Array(n);
  let minQ = 1.0;

  for (let r = n - 1; r >= 0; r--) {
    const i  = idx[r];
    // rank is r+1 (1-based); BH formula: p(r+1) * m / (r+1)
    const qi = Math.min(pvals[i] * n / (r + 1), 1.0);
    if (qi < minQ) minQ = qi;
    q[i] = minQ;
  }
  return q;
}

/** Fisher-Yates shuffle (TypedArray-compatible). */
export function shuffle(arr) {
  for (let i = arr.length - 1; i > 0; i--) {
    const j = (Math.random() * (i + 1)) | 0;
    const t = arr[i]; arr[i] = arr[j]; arr[j] = t;
  }
}

/** NES for KS: ES / mean(|same-sign null ES|). */
export function calcNES(es, nullKS) {
  const n   = nullKS.length;
  const pos = es >= 0;
  let sumSame = 0, nSame = 0, sumAll = 0;
  for (let j = 0; j < n; j++) {
    const v  = nullKS[j];
    const av = v < 0 ? -v : v;
    sumAll += av;
    if (pos ? v >= 0 : v < 0) { sumSame += av; nSame++; }
  }
  const mean = nSame > 0 ? sumSame / nSame : sumAll / n;
  return mean < 1e-12 ? 0 : es / mean;
}

/** NES for AD: AD / mean(null AD). */
export function calcNES_AD(ad, nullAD) {
  const n = nullAD.length;
  let sum = 0;
  for (let j = 0; j < n; j++) sum += nullAD[j];
  const mean = sum / n;
  return mean < 1e-12 ? 0 : ad / mean;
}

/** Two-sided empirical p for KS. Laplace-smoothed. */
export function empP_KS(obsKS, nullKS, nPerms) {
  const ab = obsKS < 0 ? -obsKS : obsKS;
  let c = 0;
  for (let j = 0; j < nPerms; j++) {
    const v = nullKS[j];
    if ((v < 0 ? -v : v) >= ab) c++;
  }
  return (c + 1) / (nPerms + 1);
}

/** One-sided empirical p for AD (larger = more extreme). */
export function empP_AD(obsAD, nullAD, nPerms) {
  let c = 0;
  for (let j = 0; j < nPerms; j++) { if (nullAD[j] >= obsAD) c++; }
  return (c + 1) / (nPerms + 1);
}

// 在 core/math.js 中添加或修改

/**
 * 核心：计算 GSEA 风格的 FDR (基于全局 NES 分布)
 * @param {number[]} obsNES - 观测到的 NES 数组 (nP)
 * @param {Float64Array[]} nullNESMat - 所有置换的 NES 矩阵 [nP][nPerms]
 * @param {boolean} twoSided - 对于 KS 需要分正负，对于 AD 只需要单边
 */
export function calcGseaFDR(obsNES, nullNESMat, twoSided = true) {
  const nP = obsNES.length;
  const nPerms = nullNESMat[0].length;
  const fdr = new Float64Array(nP);

  if (twoSided) {
    // 1. 准备全局池 (KS 逻辑：正负分开处理)
    const nullPos = [], nullNeg = [];
    for (let i = 0; i < nP; i++) {
      for (let j = 0; j < nPerms; j++) {
        const v = nullNESMat[i][j];
        if (v >= 0) nullPos.push(v); else nullNeg.push(v);
      }
    }
    // 排序以便快速计数 (可选优化)
    nullPos.sort((a, b) => a - b);
    nullNeg.sort((a, b) => a - b);

    for (let i = 0; i < nP; i++) {
      const o = obsNES[i];
      const isPos = o >= 0;
      const pool = isPos ? nullPos : nullNeg;
      const obsList = obsNES.filter(v => isPos ? v >= 0 : v < 0);

      const countNull = pool.filter(v => isPos ? v >= o : v <= o).length;
      const countObs  = obsList.filter(v => isPos ? v >= o : v <= o).length;

      const phiNull = countNull / pool.length;
      const phiObs  = countObs / obsList.length;
      fdr[i] = Math.min(1.0, phiNull / (phiObs + 1e-10));
    }
  } else {
    // 2. AD 逻辑：全是正数，单边处理
    const nullPool = [];
    for (let i = 0; i < nP; i++) {
      for (let j = 0; j < nPerms; j++) nullPool.push(nullNESMat[i][j]);
    }
    nullPool.sort((a, b) => a - b);

    for (let i = 0; i < nP; i++) {
      const o = obsNES[i];
      const countNull = nullPool.filter(v => v >= o).length;
      const countObs  = obsNES.filter(v => v >= o).length;
      
      const phiNull = countNull / nullPool.length;
      const phiObs  = countObs / nP;
      fdr[i] = Math.min(1.0, phiNull / (phiObs + 1e-10));
    }
  }
  return fdr;
}
