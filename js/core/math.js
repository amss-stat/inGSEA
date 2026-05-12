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
/**
 * GSEA-style FDR with isotonic correction.
 * Subramanian et al. 2005, Supplementary Methods.
 *
 * For positive NES threshold t:
 *   phiNull(t) = #{null NES >= t AND >= 0} / #{null NES >= 0}
 *   phiObs(t)  = #{obs  NES >= t AND >= 0} / #{obs  NES >= 0}
 *   rawFDR(t)  = phiNull(t) / phiObs(t), clamped to [0,1]
 *
 * Isotonic correction (original paper):
 *   If t' > t then FDR(t') <= FDR(t)
 *   → running minimum from least extreme to most extreme
 *
 * @param {number[]}       obsNES     observed NES [nP]
 * @param {Float64Array[]} nullNESMat null NES matrix [nP][nPerms]
 * @param {boolean}        twoSided   true=KS(pos/neg split), false=AD(all pos)
 * @returns {Float64Array} fdr [nP], isotonically corrected
 */
export function calcGseaFDR(obsNES, nullNESMat, twoSided) {
  const nP     = obsNES.length;
  const nPerms = nullNESMat[0].length;
  const fdr    = new Float64Array(nP).fill(1.0);

  // ── Binary search helpers (array sorted ascending) ──────────
  // Count elements >= threshold
  function countGeq(sorted, t) {
    let lo = 0, hi = sorted.length;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (sorted[mid] < t) lo = mid + 1; else hi = mid;
    }
    return sorted.length - lo;
  }

  // Count elements <= threshold
  function countLeq(sorted, t) {
    let lo = 0, hi = sorted.length;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (sorted[mid] <= t) lo = mid + 1; else hi = mid;
    }
    return lo;
  }

  /**
   * Process one side (positive or negative).
   *
   * @param {number[]} idxSide  indices into obsNES for this side
   * @param {number[]} nullPool sorted ascending, null NES for this side
   * @param {boolean}  geq      true=positive side, false=negative side
   */
  function processSide(idxSide, nullPool, geq) {
    const nObs  = idxSide.length;
    const nNull = nullPool.length;
    if (nObs === 0 || nNull === 0) return;

    // Sort: most extreme first
    // Positive: descending (largest NES first, r=0)
    // Negative: ascending  (most negative first, r=0)
    idxSide.sort((a, b) =>
      geq ? obsNES[b] - obsNES[a]
          : obsNES[a] - obsNES[b]
    );

    // ── Step 1: compute rawFDR for each rank ─────────────────
    // At rank r (0-based, most extreme first):
    //   phiObs(r)  = (r+1) / nObs
    //   phiNull(r) = countExtreme(nullPool, t) / nNull
    const rawFDR = new Float64Array(nObs);

    for (let r = 0; r < nObs; r++) {
      const t = obsNES[idxSide[r]];

      const nullCount = geq
        ? countGeq(nullPool, t)
        : countLeq(nullPool, t);

      const phiNull = nullCount / nNull;
      const phiObs  = (r + 1)  / nObs;

      rawFDR[r] = Math.min(phiNull / phiObs, 1.0);
    }

    // ── Step 2: isotonic correction ───────────────────────────
    // Original paper: "if t' > t then FDR(t') <= FDR(t)"
    // In our ordering: r=0 is most extreme (t' > t for smaller r)
    // So FDR must be non-increasing as r decreases.
    // Equivalently: FDR is non-decreasing as r increases.
    //
    // Correct direction: scan from r=nObs-1 DOWN to r=0,
    // carrying running minimum rightward→leftward.
    //
    //   corrFDR[nObs-1] = rawFDR[nObs-1]
    //   corrFDR[r]      = min(rawFDR[r], corrFDR[r+1])
    //
    // This ensures corrFDR[r] <= corrFDR[r+1] for all r,
    // i.e. more extreme NES gets <= FDR than less extreme NES.

    let minFDR = 1.0;
    for (let r = nObs - 1; r >= 0; r--) {
      if (rawFDR[r] < minFDR) minFDR = rawFDR[r];
      fdr[idxSide[r]] = minFDR;
    }
  }

  if (twoSided) {
    // ── Build null pools ──────────────────────────────────────
    const nullPos = [];
    const nullNeg = [];
    for (let pi = 0; pi < nP; pi++) {
      const vec = nullNESMat[pi];
      for (let j = 0; j < nPerms; j++) {
        const v = vec[j];
        if (v >= 0) nullPos.push(v);
        else        nullNeg.push(v);
      }
    }
    nullPos.sort((a, b) => a - b);
    nullNeg.sort((a, b) => a - b);

    // ── Split observed by sign ────────────────────────────────
    const idxPos = [], idxNeg = [];
    for (let i = 0; i < nP; i++) {
      if (obsNES[i] >= 0) idxPos.push(i);
      else                idxNeg.push(i);
    }

    processSide(idxPos, nullPos, true);
    processSide(idxNeg, nullNeg, false);

  } else {
    // ── Single-sided (AD: all positive) ──────────────────────
    const nullPool = [];
    for (let pi = 0; pi < nP; pi++) {
      const vec = nullNESMat[pi];
      for (let j = 0; j < nPerms; j++) nullPool.push(vec[j]);
    }
    nullPool.sort((a, b) => a - b);

    const idxAll = Array.from({ length: nP }, (_, i) => i);
    processSide(idxAll, nullPool, true);
  }

  return fdr;
}


/**
 * Permutation-based FDR on pCauchy.
 *
 * Uses the null distribution of pCauchy values derived directly
 * from permutation ES and AD statistics, preserving the
 * inter-pathway correlation structure.
 *
 * Algorithm:
 *   1. For each permutation j, compute null pCauchy(pi,j) for
 *      all pathways using rank-based p-value approximation.
 *   2. For each observed pCauchy threshold t:
 *        V(t) = expected # pathways with null pCauchy <= t
 *             = (1/nPerms) * Σ_j #{pi : null_pCauchy(pi,j) <= t}
 *        R(t) = # pathways with obs pCauchy <= t
 *        FDR(t) = V(t) / R(t), clamped to [0,1]
 *   3. Isotonic correction: running minimum from smallest to
 *      largest pCauchy (most to least significant).
 *
 * Correlation handling:
 *   Because all pathways share the same permutation j,
 *   the joint null distribution reflects actual inter-pathway
 *   correlations without any independence assumption.
 *
 * Self-comparison bias:
 *   Rank-based null p-values include the observation itself
 *   in the reference distribution, causing slight conservative
 *   bias of order 1/nPerms. Acceptable for nPerms >= 500.
 *
 * @param {number[]}       obsPCauchy   observed pCauchy [nP]
 * @param {Float64Array[]} nullKS       null KS ES matrix [nP][nPerms]
 * @param {Float64Array[]} nullAD       null AD matrix    [nP][nPerms]
 * @param {number}         nPerms
 * @returns {Float64Array} fdr [nP], isotonically corrected
 */
export function permFDR(obsPCauchy, nullKS, nullAD, nPerms) {
  const nP  = obsPCauchy.length;
  const fdr = new Float64Array(nP).fill(1.0);
  if (nP < 2) return fdr;

  // ── Step 1: Pre-sort null KS and null AD for rank lookup ───
  // For each pathway pi, we need the rank of nullKS[pi][j]
  // within nullKS[pi] to get the null p-value.
  //
  // Approach: sort a copy of each row, then use binary search
  // to find the rank of each element.

  // sortedKS[pi] = nullKS[pi] sorted ascending (for |KS|, two-sided)
  // sortedAD[pi] = nullAD[pi] sorted ascending (for AD, one-sided)
  const sortedAbsKS = new Array(nP);
  const sortedAD    = new Array(nP);

  for (let pi = 0; pi < nP; pi++) {
    // KS is two-sided: use |KS|
    const absKS = new Float64Array(nPerms);
    const ksVec = nullKS[pi];
    for (let j = 0; j < nPerms; j++)
      absKS[j] = ksVec[j] < 0 ? -ksVec[j] : ksVec[j];
    absKS.sort();
    sortedAbsKS[pi] = absKS;

    // AD is one-sided (always positive)
    const adCopy = nullAD[pi].slice().sort();
    sortedAD[pi] = adCopy;
  }

  // ── Step 2: Build null pCauchy matrix ──────────────────────
  // nullPCauchy[j] = Float64Array of length nP
  // (one value per pathway for permutation j)
  //
  // We don't store the full matrix to save memory.
  // Instead, for each threshold t we need to count
  // #{pi : nullPCauchy(pi,j) <= t} across all j.
  //
  // Strategy: build the full null pCauchy distribution as a
  // flat sorted array of nP*nPerms values, then use binary
  // search for FDR counting.
  //
  // Memory: nP * nPerms * 8 bytes = 500 * 1000 * 8 = 4MB — OK

  // Helper: rank-based p-value (Laplace smoothed)
  // countGeqAbs: # values in sortedAbs >= |obs| (two-sided)
  function nullPKS(pi, absObs) {
    // Binary search in sortedAbsKS[pi] (ascending)
    // count >= absObs
    const arr = sortedAbsKS[pi];
    let lo = 0, hi = nPerms;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (arr[mid] < absObs) lo = mid + 1; else hi = mid;
    }
    const c = nPerms - lo;
    return (c + 1) / (nPerms + 1);
  }

  // countGeq: # values in sortedAD >= obs (one-sided)
  function nullPAD(pi, obs) {
    const arr = sortedAD[pi];
    let lo = 0, hi = nPerms;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (arr[mid] < obs) lo = mid + 1; else hi = mid;
    }
    const c = nPerms - lo;
    return (c + 1) / (nPerms + 1);
  }

  // Build flat null pCauchy array
  const nullPCauchyFlat = new Float64Array(nP * nPerms);
  let ptr = 0;

  for (let j = 0; j < nPerms; j++) {
    for (let pi = 0; pi < nP; pi++) {
      const ksVal = nullKS[pi][j];
      const adVal = nullAD[pi][j];

      const absKS = ksVal < 0 ? -ksVal : ksVal;
      const pKS   = nullPKS(pi, absKS);
      const pAD   = nullPAD(pi, adVal);

      nullPCauchyFlat[ptr++] = cauchyCombine(pKS, pAD);
    }
  }

  // Sort ascending for binary search
  nullPCauchyFlat.sort();
  const nNull = nullPCauchyFlat.length;  // = nP * nPerms

  // ── Step 3: FDR for each observed pCauchy ──────────────────
  // Sort observed pCauchy ascending (smallest = most significant first)
  const idxObs = Array.from({ length: nP }, (_, i) => i);
  idxObs.sort((a, b) => obsPCauchy[a] - obsPCauchy[b]);

  // For threshold t = obsPCauchy[idxObs[r]]:
  //   R(t) = r + 1   (# obs with pCauchy <= t, since sorted ascending)
  //   V(t) = #{null pCauchy <= t} / nPerms
  //          (expected # pathways per permutation with null pCauchy <= t)
  //
  // FDR(t) = V(t) / R(t)

  // Count null pCauchy <= t via binary search (upper_bound)
  function countNullLeq(t) {
    let lo = 0, hi = nNull;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (nullPCauchyFlat[mid] <= t) lo = mid + 1; else hi = mid;
    }
    return lo;  // # elements <= t
  }

  const rawFDR = new Float64Array(nP);
  for (let r = 0; r < nP; r++) {
    const idx = idxObs[r];
    const t   = obsPCauchy[idx];

    // V(t): expected # null pCauchy <= t per permutation
    // = (total null pCauchy <= t) / nPerms
    const nullLeq = countNullLeq(t);
    const Vt = nullLeq / nPerms;  // expected false positives

    // R(t): observed # pCauchy <= t
    const Rt = r + 1;

    rawFDR[r] = Math.min(Vt / Rt, 1.0);
  }

  // ── Step 4: Isotonic correction ────────────────────────────
  // FDR must be non-decreasing as pCauchy increases
  // (less significant thresholds must have >= FDR)
  //
  // We've sorted ascending pCauchy: r=0 most significant.
  // Correction: scan from r=nP-1 (least significant) down to r=0,
  // carrying running minimum — ensures FDR non-decreasing with p.
  //
  // After correction: fdr[r=0] <= fdr[r=1] <= ... <= fdr[r=nP-1]

  let minFDR = 1.0;
  for (let r = nP - 1; r >= 0; r--) {
    if (rawFDR[r] < minFDR) minFDR = rawFDR[r];
    fdr[idxObs[r]] = minFDR;
  }

  return fdr;
}


/**
 * Storey (2002) q-value.
 *
 * Estimates π₀ (proportion of true nulls) from the p-value
 * distribution, then applies a π₀-weighted BH correction.
 *
 * π₀ estimation strategy:
 *   π₀(λ) = #{p > λ} / (m × (1 − λ))
 *   Evaluate at λ ∈ {0.05, 0.10, …, 0.50}.
 *   Use the estimate at the largest λ where at least
 *   max(5, m×0.10) p-values exceed λ (stability criterion).
 *   If no stable λ found: π₀ = 1 (reduces to BH).
 *
 * For m < 10: always uses π₀ = 1 (too few p-values to
 * estimate π₀ reliably).
 *
 * Monotonicity: q-values are guaranteed non-decreasing with
 * p-values via right-to-left running minimum (same as BH).
 *
 * @param  {number[]}      pvals  raw p-values (any order)
 * @returns {Float64Array}        q-values aligned to input
 */
export function storeyQ(pvals) {
  const m = pvals.length;
  if (m === 0) return new Float64Array(0);
  if (m === 1) return Float64Array.from([Math.min(pvals[0], 1)]);

  // ── Step 1: estimate π₀ ──────────────────────────────────
  let pi0 = 1.0;

  if (m >= 10) {
    const lambdas = [0.05, 0.10, 0.15, 0.20, 0.25,
                     0.30, 0.35, 0.40, 0.45, 0.50];
    const minAbove = Math.max(5, Math.floor(m * 0.10));

    // Scan from largest λ downward; take first stable estimate
    for (let li = lambdas.length - 1; li >= 0; li--) {
      const lam   = lambdas[li];
      const above = pvals.filter(p => p > lam).length;
      if (above < minAbove) continue;
      pi0 = Math.min(above / (m * (1 - lam)), 1.0);
      break;
    }
    // If no stable λ found, pi0 remains 1.0 (conservative = BH)
  }

  // ── Step 2: π₀-weighted BH ───────────────────────────────
  // q(i) = min_{j ≥ i} { π₀ × m × p(j) / j }   (1-based rank j)
  // Right-to-left running minimum ensures monotonicity.
  const idx = Array.from({ length: m }, (_, i) => i);
  idx.sort((a, b) => pvals[a] - pvals[b]);   // ascending p-value

  const q    = new Float64Array(m);
  let   minQ = 1.0;

  for (let r = m - 1; r >= 0; r--) {
    const i  = idx[r];
    const qi = Math.min(pi0 * pvals[i] * m / (r + 1), 1.0);
    if (qi < minQ) minQ = qi;
    q[i] = minQ;
  }

  return q;
}
