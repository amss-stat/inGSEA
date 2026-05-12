// ═══════════════════════════════════════════════════════════
//  core/math.js  ·  Linear-algebra primitives & statistics
// ═══════════════════════════════════════════════════════════
'use strict';

/**
 * Welch SNR (signal-to-noise): (μ1−μ2)/(σ1+σ2)
 * Returns Float64Array[nGenes].
 */
export function calcSNR(mat, nGenes, cIdx, tIdx) {
  const nc = cIdx.length, nt = tIdx.length;
  const snr = new Float64Array(nGenes);
  for (let g = 0; g < nGenes; g++) {
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
 * Descending argsort of snr array.
 * Returns Int32Array[nGenes].
 */
export function rankOrder(snr, nGenes) {
  const ord = new Int32Array(nGenes);
  for (let i = 0; i < nGenes; i++) ord[i] = i;
  ord.sort((a, b) => snr[b] - snr[a]);
  return ord;
}

/**
 * GSEA enrichment walk (KS-style ES + Anderson–Darling statistic).
 *
 * @param {Float64Array} snr   – signed SNR, length nGenes
 * @param {Int32Array}   ord   – gene indices sorted by descending SNR
 * @param {Uint8Array}   mask  – 1 if gene is in set, else 0
 * @param {number}       nGenes
 * @param {boolean}      wantCurve  – if true, allocate and return curve
 * @returns {{ ks:number, ad:number, curve:Float64Array|null, peakIdx:number }}
 */
export function calcGSEAStats(snr, ord, mask, nGenes, wantCurve) {
  let nHits = 0, hitSum = 0;
  for (let i = 0; i < nGenes; i++) {
    if (mask[ord[i]]) { nHits++; hitSum += Math.abs(snr[ord[i]]); }
  }
  if (nHits === 0 || nHits === nGenes) {
    return {
      ks: 0, ad: 0, peakIdx: 0,
      curve: wantCurve ? new Float64Array(nGenes) : null
    };
  }

  const nMiss    = nGenes - nHits;
  const missStep = -1 / nMiss;
  const invHit   = 1 / (hitSum + 1e-12);
  const curve    = wantCurve ? new Float64Array(nGenes) : null;

  let cum = 0, ks = 0, maxAbs = 0, peakIdx = 0, ad = 0;

  for (let i = 0; i < nGenes; i++) {
    cum += mask[ord[i]] ? Math.abs(snr[ord[i]]) * invHit : missStep;
    if (wantCurve) curve[i] = cum;

    const a = Math.abs(cum);
    if (a > maxAbs) { maxAbs = a; ks = cum; peakIdx = i; }

    // Anderson–Darling: Σ ES²/(F·(1-F)) for i = 0..nGenes-2
    if (i < nGenes - 1) {
      const F = (i + 1) / nGenes;
      ad += (cum * cum) / (F * (1 - F));
    }
  }

  return { ks, ad, curve, peakIdx };
}

/**
 * Cauchy combination of p-values (Liu & Xie 2020).
 * Handles p=0 and p=1 edge cases.
 */
export function cauchyCombine(...pvals) {
  let s = 0;
  for (const p of pvals) {
    const pc = Math.max(Math.min(p, 1 - 1e-16), 1e-16);
    s += Math.tan((0.5 - pc) * Math.PI);
  }
  s /= pvals.length;
  return 0.5 - Math.atan(s) / Math.PI;
}

/**
 * Benjamini–Hochberg FDR from sorted (ascending) p-values array.
 * Returns Float64Array of q-values aligned to original order.
 */
export function bhFDR(pvals) {
  const n = pvals.length;
  const idx = Array.from({ length: n }, (_, i) => i);
  idx.sort((a, b) => pvals[a] - pvals[b]);
  const q = new Float64Array(n);
  let minQ = 1;
  for (let r = n - 1; r >= 0; r--) {
    const i = idx[r];
    minQ = Math.min(pvals[i] * n / (r + 1), minQ);
    q[i] = minQ;
  }
  return q;
}

/**
 * Fisher–Yates shuffle in place.
 */
export function shuffle(arr) {
  for (let i = arr.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    const t = arr[i]; arr[i] = arr[j]; arr[j] = t;
  }
}

/**
 * Normalised Enrichment Score.
 * NES = ES / mean(|null_ES with same sign|)
 * Falls back to ES / mean(|all null_ES|) if sign pool is empty.
 */
export function calcNES(es, nullES) {
  const sameSign = nullES.filter(v => (es >= 0 ? v >= 0 : v < 0));
  const pool     = sameSign.length > 0 ? sameSign : nullES;
  const mean     = pool.reduce((s, v) => s + Math.abs(v), 0) / (pool.length || 1);
  return mean < 1e-12 ? 0 : es / mean;
}

/**
 * NES for AD (use absolute null AD since AD ≥ 0).
 */
export function calcNES_AD(ad, nullAD) {
  const mean = nullAD.reduce((s, v) => s + v, 0) / (nullAD.length || 1);
  return mean < 1e-12 ? 0 : ad / mean;
}
