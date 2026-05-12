// ═══════════════════════════════════════════════════════════
//  core/math.js  ·  Pure numeric primitives
//  Performance: calcSNR accepts pre-allocated output buffer
// ═══════════════════════════════════════════════════════════
'use strict';

/**
 * Welch SNR: (μ_case − μ_ctrl) / (σ_case + σ_ctrl + ε)
 *
 * @param {Float64Array[]} mat
 * @param {number}         nG
 * @param {ArrayLike}      cIdx  case indices (Uint16Array view OK)
 * @param {ArrayLike}      tIdx  ctrl indices
 * @param {Float64Array}   [out] pre-allocated output buffer (length nG)
 * @returns {Float64Array}  out (or new array if not provided)
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
 * Descending argsort into pre-allocated Int32Array.
 * @param {Float64Array} snr
 * @param {number}       nG
 * @param {Int32Array}   ord  pre-allocated, length nG
 */
export function rankOrder(snr, nG, ord) {
  for (let i = 0; i < nG; i++) ord[i] = i;
  ord.sort((a, b) => snr[b] - snr[a]);
}

/**
 * GSEA enrichment walk.
 * Returns { ks, ad, curve (or null), peakIdx }.
 */
export function calcGSEAStats(snr, ord, mask, nG, wantCurve) {
  let nHits = 0, hitSum = 0;
  for (let i = 0; i < nG; i++) {
    if (mask[ord[i]]) {
      nHits++;
      hitSum += snr[ord[i]] < 0 ? -snr[ord[i]] : snr[ord[i]];
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
    const absSnr = snr[ord[i]];
    cum += mask[ord[i]]
      ? (absSnr < 0 ? -absSnr : absSnr) * invHit
      : missStep;

    if (wantCurve) curve[i] = cum;

    const a = cum < 0 ? -cum : cum;
    if (a > maxAbs) { maxAbs = a; ks = cum; peakIdx = i; }

    if (i < nG - 1) {
      const F = (i + 1) / nG;
      ad += (cum * cum) / (F * (1.0 - F));
    }
  }

  return { ks, ad, curve, peakIdx };
}

/** Cauchy combination of p-values. */
export function cauchyCombine(...pvals) {
  let s = 0;
  for (let p of pvals) {
    p = p < 1e-16 ? 1e-16 : p > 1 - 1e-16 ? 1 - 1e-16 : p;
    s += Math.tan((0.5 - p) * Math.PI);
  }
  s /= pvals.length;
  return 0.5 - Math.atan(s) / Math.PI;
}

/** Benjamini-Hochberg FDR. */
export function bhFDR(pvals) {
  const n   = pvals.length;
  const idx = Array.from({ length: n }, (_, i) => i);
  idx.sort((a, b) => pvals[a] - pvals[b]);
  const q = new Float64Array(n);
  let min = 1.0;
  for (let r = n - 1; r >= 0; r--) {
    const i  = idx[r];
    const qi = pvals[i] * n / (r + 1);
    min   = qi < min ? qi : min;
    q[i]  = min;
  }
  return q;
}

/** Fisher-Yates shuffle (works for TypedArray). */
export function shuffle(arr) {
  for (let i = arr.length - 1; i > 0; i--) {
    const j = (Math.random() * (i + 1)) | 0;
    const t = arr[i]; arr[i] = arr[j]; arr[j] = t;
  }
}

/** NES for KS statistic. */
export function calcNES(es, nullKS) {
  const n   = nullKS.length;
  const pos = es >= 0;
  let sumSame = 0, nSame = 0, sumAll = 0;
  for (let j = 0; j < n; j++) {
    const v = nullKS[j];
    const av = v < 0 ? -v : v;
    sumAll += av;
    if (pos ? v >= 0 : v < 0) { sumSame += av; nSame++; }
  }
  const mean = nSame > 0 ? sumSame / nSame : sumAll / n;
  return mean < 1e-12 ? 0 : es / mean;
}

/** NES for AD statistic. */
export function calcNES_AD(ad, nullAD) {
  const n = nullAD.length;
  let sum = 0;
  for (let j = 0; j < n; j++) sum += nullAD[j];
  const mean = sum / n;
  return mean < 1e-12 ? 0 : ad / mean;
}

/** Empirical p for KS (two-sided). */
export function empP_KS(obsKS, nullKS, nPerms) {
  const ab = obsKS < 0 ? -obsKS : obsKS;
  let c = 0;
  for (let j = 0; j < nPerms; j++) {
    const v = nullKS[j];
    if ((v < 0 ? -v : v) >= ab) c++;
  }
  return (c + 1) / (nPerms + 1);
}

/** Empirical p for AD (one-sided). */
export function empP_AD(obsAD, nullAD, nPerms) {
  let c = 0;
  for (let j = 0; j < nPerms; j++) { if (nullAD[j] >= obsAD) c++; }
  return (c + 1) / (nPerms + 1);
}
