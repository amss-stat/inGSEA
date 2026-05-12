// ═══════════════════════════════════════════════════════════
//  core/math.js  ·  Pure numeric primitives
//  NO distribution fitting here — that lives in webr-bridge.js
// ═══════════════════════════════════════════════════════════
'use strict';

// ── Signal-to-noise ratio ────────────────────────────────────
/**
 * Welch SNR: (μ_case − μ_ctrl) / (σ_case + σ_ctrl)
 *
 * cIdx / tIdx may be plain arrays OR Uint16Array views —
 * we iterate with numeric indexing, which works for both.
 *
 * @param {Float64Array[]} mat   genes × samples
 * @param {number}         nG   number of genes
 * @param {ArrayLike}      cIdx case sample indices
 * @param {ArrayLike}      tIdx ctrl sample indices
 * @returns {Float64Array}
 */
export function calcSNR(mat, nG, cIdx, tIdx) {
  const nc = cIdx.length, nt = tIdx.length;
  const snr = new Float64Array(nG);
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

// ── Descending argsort ───────────────────────────────────────
/**
 * Return Int32Array of indices sorted by descending snr[i].
 * @param {Float64Array} snr
 * @param {number}       nG
 * @returns {Int32Array}
 */
export function rankOrder(snr, nG) {
  const ord = new Int32Array(nG);
  for (let i = 0; i < nG; i++) ord[i] = i;
  // Standard sort on typed array — V8 uses TimSort, stable
  ord.sort((a, b) => snr[b] - snr[a]);
  return ord;
}

// ── GSEA enrichment walk ─────────────────────────────────────
/**
 * Compute KS-style ES and Anderson–Darling accumulation.
 *
 * ES walk:
 *   • hit gene  → +|snr[gene]| / hitSum
 *   • miss gene → −1 / nMiss
 * KS  = running max |cum| with sign
 * AD  = Σ_{i=0}^{n-2}  cum² / (F_i · (1 − F_i))
 *
 * @param {Float64Array} snr
 * @param {Int32Array}   ord      descending SNR order
 * @param {Uint8Array}   mask     1 = gene in set
 * @param {number}       nG
 * @param {boolean}      wantCurve  allocate + return curve array
 * @returns {{ ks, ad, curve, peakIdx }}
 */
export function calcGSEAStats(snr, ord, mask, nG, wantCurve) {
  // Pre-compute hit count and weighted sum
  let nHits = 0, hitSum = 0;
  for (let i = 0; i < nG; i++) {
    if (mask[ord[i]]) {
      nHits++;
      hitSum += Math.abs(snr[ord[i]]);
    }
  }

  // Degenerate cases → ES = 0
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
    // Accumulate running enrichment score
    cum += mask[ord[i]] ? Math.abs(snr[ord[i]]) * invHit : missStep;
    if (wantCurve) curve[i] = cum;

    // KS: track signed extremum
    const a = cum < 0 ? -cum : cum;   // Math.abs without function call
    if (a > maxAbs) { maxAbs = a; ks = cum; peakIdx = i; }

    // AD: use interior points only (i = 0 … nG-2)
    if (i < nG - 1) {
      const F = (i + 1) / nG;
      ad += (cum * cum) / (F * (1.0 - F));
    }
  }

  return { ks, ad, curve, peakIdx };
}

// ── Cauchy combination ───────────────────────────────────────
/**
 * Combine an arbitrary number of p-values via Cauchy weights.
 * Liu & Xie (2020), equal weights.
 * Handles p ≈ 0 and p ≈ 1 safely.
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

// ── Benjamini–Hochberg FDR ───────────────────────────────────
/**
 * BH-adjusted q-values for an array of p-values (any order).
 * Returns Float64Array aligned to input order.
 *
 * @param {number[]} pvals
 * @returns {Float64Array}
 */
export function bhFDR(pvals) {
  const n   = pvals.length;
  const idx = Array.from({ length: n }, (_, i) => i);
  idx.sort((a, b) => pvals[a] - pvals[b]);   // sort indices by p ascending

  const q = new Float64Array(n);
  let min = 1.0;
  // Sweep from largest rank downward
  for (let r = n - 1; r >= 0; r--) {
    const i  = idx[r];
    const qi = pvals[i] * n / (r + 1);
    min    = qi < min ? qi : min;
    q[i]   = min;
  }
  return q;
}

// ── Fisher–Yates shuffle ─────────────────────────────────────
/**
 * In-place shuffle of any array-like with numeric indices.
 * Works for both plain Array and TypedArray.
 */
export function shuffle(arr) {
  for (let i = arr.length - 1; i > 0; i--) {
    const j = (Math.random() * (i + 1)) | 0;   // bitwise OR = Math.floor for positive
    const t = arr[i]; arr[i] = arr[j]; arr[j] = t;
  }
}

// ── Normalised Enrichment Scores ─────────────────────────────
/**
 * NES for KS (standard GSEA normalisation).
 * Divide observed ES by mean of same-sign null ES values.
 * Falls back to mean of |all| null ES if same-sign pool is empty.
 *
 * Uses typed array iteration throughout — no filter/map.
 *
 * @param {number}       es       observed ES
 * @param {Float64Array} nullKS   null distribution of KS statistics
 * @returns {number}
 */
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

/**
 * NES for AD.
 * AD is always ≥ 0, so normalise by mean of all null AD values.
 *
 * @param {number}       ad       observed AD
 * @param {Float64Array} nullAD   null distribution of AD statistics
 * @returns {number}
 */
export function calcNES_AD(ad, nullAD) {
  const n = nullAD.length;
  let sum = 0;
  for (let j = 0; j < n; j++) sum += nullAD[j];
  const mean = sum / n;
  return mean < 1e-12 ? 0 : ad / mean;
}

// ── Empirical p-values ───────────────────────────────────────
/**
 * Two-sided empirical p for KS: P(|null| >= |obs|).
 * Uses Laplace smoothing: (count + 1) / (nPerms + 1).
 *
 * @param {number}       obsKS
 * @param {Float64Array} nullKS
 * @param {number}       nPerms
 */
export function empP_KS(obsKS, nullKS, nPerms) {
  const ab = obsKS < 0 ? -obsKS : obsKS;
  let c = 0;
  for (let j = 0; j < nPerms; j++) {
    const v = nullKS[j]; if ((v < 0 ? -v : v) >= ab) c++;
  }
  return (c + 1) / (nPerms + 1);
}

/**
 * One-sided empirical p for AD: P(null >= obs).
 *
 * @param {number}       obsAD
 * @param {Float64Array} nullAD
 * @param {number}       nPerms
 */
export function empP_AD(obsAD, nullAD, nPerms) {
  let c = 0;
  for (let j = 0; j < nPerms; j++) { if (nullAD[j] >= obsAD) c++; }
  return (c + 1) / (nPerms + 1);
}
