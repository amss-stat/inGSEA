// ═══════════════════════════════════════════════════════════
//  core/distributions.js
//
//  Parametric null-distribution fitting using jStat for
//  the regularised incomplete gamma function (the only
//  special function we need that is hard to implement
//  correctly from scratch).
//
//  jStat.gammainc(x, a) = P(a, x)  [lower regularised]
//  jStat.lngamma(a)     = log Γ(a)
//
//  GG parameterisation: (μ, σ, Q) — identical to flexsurv.
//  Gamma parameterisation: (shape k, rate r).
//
//  Optimizer: Nelder-Mead simplex (pure JS, ~100 lines).
//  This is the same algorithm R's optim(method="Nelder-Mead")
//  uses. It does not require gradients and is robust for
//  low-dimensional MLE problems (2–3 parameters).
// ═══════════════════════════════════════════════════════════
'use strict';

// ── jStat availability check ─────────────────────────────────
function _jstat() {
  if (typeof jStat === 'undefined')
    throw new Error('jStat not loaded — parametric engine unavailable');
  return jStat;
}

// ── Nelder-Mead optimizer ────────────────────────────────────
/**
 * Minimise f(x) starting from x0 using the Nelder-Mead simplex.
 * Pure JS, no dependencies.
 *
 * @param {(x: number[]) => number} f
 * @param {number[]} x0   initial parameter vector
 * @param {object}   opts
 * @returns {{ x: number[], fval: number, converged: boolean }}
 */
export function nelderMead(f, x0, opts = {}) {
  const {
    maxIter = 2000,
    tol     = 1e-10,
    alpha   = 1.0,   // reflection
    gamma   = 2.0,   // expansion
    rho     = 0.5,   // contraction
    sigma   = 0.5    // shrink
  } = opts;

  const n = x0.length;

  // Build initial simplex: x0 plus n perturbed vertices
  let simplex = [x0.slice()];
  for (let i = 0; i < n; i++) {
    const v = x0.slice();
    v[i] += (Math.abs(v[i]) > 1e-8) ? 0.05 * v[i] : 0.00025;
    simplex.push(v);
  }

  let fvals = simplex.map(v => f(v));

  for (let iter = 0; iter < maxIter; iter++) {
    // Sort by function value
    const idx = Array.from({ length: n + 1 }, (_, i) => i);
    idx.sort((a, b) => fvals[a] - fvals[b]);
    simplex = idx.map(i => simplex[i]);
    fvals   = idx.map(i => fvals[i]);

    // Convergence check
    if (fvals[n] - fvals[0] < tol) {
      return { x: simplex[0], fval: fvals[0], converged: true };
    }

    // Centroid of all vertices except worst
    const cent = new Array(n).fill(0);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        cent[j] += simplex[i][j] / n;

    // Reflect worst through centroid
    const xr = cent.map((c, j) => c + alpha * (c - simplex[n][j]));
    const fr = f(xr);

    if (fr < fvals[0]) {
      // Try expansion
      const xe = cent.map((c, j) => c + gamma * (xr[j] - c));
      const fe = f(xe);
      if (fe < fr) { simplex[n] = xe; fvals[n] = fe; }
      else         { simplex[n] = xr; fvals[n] = fr; }
    } else if (fr < fvals[n - 1]) {
      simplex[n] = xr; fvals[n] = fr;
    } else {
      // Contraction
      const xc = cent.map((c, j) => c + rho * (simplex[n][j] - c));
      const fc = f(xc);
      if (fc < fvals[n]) {
        simplex[n] = xc; fvals[n] = fc;
      } else {
        // Shrink all vertices toward best
        for (let i = 1; i <= n; i++) {
          simplex[i] = simplex[i].map((v, j) =>
            simplex[0][j] + sigma * (v - simplex[0][j]));
          fvals[i] = f(simplex[i]);
        }
      }
    }
  }

  return { x: simplex[0], fval: fvals[0], converged: false };
}

// ── Generalised Gamma distribution ───────────────────────────
//
// Parameterisation (flexsurv / Prentice 1974):
//   w     = (log t - μ) / σ
//   k     = 1 / Q²          (shape)
//   u     = k · exp(Q · w)
//
// log f(t; μ,σ,Q) =
//   log|Q| + k·log k − log Γ(k) − log σ − log t + k·Q·w − u
//
// Special case Q→0: log-Normal(μ, σ)

/**
 * GG log-pdf (scalar, returns -Infinity for invalid inputs).
 */
function ggLogPdf(t, mu, sigma, Q) {
  if (t <= 0 || sigma <= 0) return -Infinity;
  const js = _jstat();

  const lnt = Math.log(t);

  if (Math.abs(Q) < 1e-8) {
    // Log-normal limit
    const z = (lnt - mu) / sigma;
    return -lnt - Math.log(sigma) - 0.5 * Math.log(2 * Math.PI) - 0.5 * z * z;
  }

  const k = 1.0 / (Q * Q);
  const w = (lnt - mu) / sigma;
  const u = k * Math.exp(Q * w);

  // Catch numerical overflow in u
  if (!isFinite(u) || u < 0) return -Infinity;

  return Math.log(Math.abs(Q))
    + k * Math.log(k)
    - js.lngamma(k)
    - Math.log(sigma)
    - lnt
    + k * Q * w
    - u;
}

/**
 * GG negative log-likelihood over a data array.
 * Parameters: par = [mu, log(sigma), Q]  (log-transform ensures σ > 0)
 */
function ggNegLL(par, data) {
  const mu    = par[0];
  const sigma = Math.exp(par[1]);
  const Q     = par[2];
  if (sigma < 1e-8 || sigma > 1e8) return 1e20;
  let ll = 0;
  for (let i = 0; i < data.length; i++) {
    const v = ggLogPdf(data[i], mu, sigma, Q);
    if (!isFinite(v)) return 1e20;
    ll += v;
  }
  return -ll;
}

/**
 * GG survival: P(T > x | μ, σ, Q)
 * Uses jStat.gammainc for the regularised incomplete gamma.
 *
 * jStat.gammainc(x, a) = P(a, x) = lower regularised gamma
 * So upper tail = 1 - jStat.gammainc(u, k)
 */
function ggSurvival(x, mu, sigma, Q) {
  if (x <= 0) return 1;
  const js = _jstat();

  if (Math.abs(Q) < 1e-8) {
    // Log-normal: S(x) = Φ(-(log(x)-μ)/σ) = 1 - Φ(z)
    const z = (Math.log(x) - mu) / sigma;
    return 1 - js.normal.cdf(z, 0, 1);
  }

  const k = 1.0 / (Q * Q);
  const w = (Math.log(x) - mu) / sigma;
  const u = k * Math.exp(Q * w);

  if (!isFinite(u) || u < 0) return Q > 0 ? 1 : 0;

  // jStat.gammainc(u, k) = P(k, u) = lower regularised gamma
  const lowerP = js.gammainc(u, k);

  if (Q > 0) {
    // S(x) = P(Gamma(k,1) > u) = 1 - P(k,u)
    return Math.max(0, Math.min(1, 1 - lowerP));
  } else {
    // Q < 0: direction reverses
    // S(x) = P(Gamma(k,1) < u) = P(k,u)
    return Math.max(0, Math.min(1, lowerP));
  }
}

/**
 * Fit GG to positive data using Nelder-Mead MLE.
 *
 * @param {Float64Array | number[]} data  positive values
 * @returns {{ mu, sigma, Q, ok: boolean }}
 */
export function fitGenGamma(data) {
  // Filter positive values only
  const pos = [];
  for (let i = 0; i < data.length; i++) {
    if (data[i] > 0) pos.push(data[i]);
  }
  if (pos.length < 5) return { ok: false };

  // Seed from log-normal moments (Q=0 starting point)
  const logD = pos.map(v => Math.log(v));
  const muS  = logD.reduce((s, v) => s + v, 0) / pos.length;
  const varS = logD.reduce((s, v) => s + (v - muS) ** 2, 0) / (pos.length - 1);
  const sigS = Math.max(Math.sqrt(varS), 1e-4);

  const result = nelderMead(
    par => ggNegLL(par, pos),
    [muS, Math.log(sigS), 0.1],
    { maxIter: 3000, tol: 1e-10 }
  );

  const mu    = result.x[0];
  const sigma = Math.exp(result.x[1]);
  const Q     = result.x[2];

  if (!isFinite(result.fval) || sigma < 1e-8) return { ok: false };

  return { mu, sigma, Q, ok: true };
}

/**
 * Compute upper-tail p-value from fitted GG model.
 * Scales observation by same sFac used in fitting.
 *
 * @param {number}  obsAD      observed AD statistic
 * @param {Float64Array} nullAD  null AD distribution (unscaled)
 * @param {number}  empP       empirical p-value (fallback)
 * @returns {number}
 */
export function pvalGG(obsAD, nullAD, empP) {
  try {
    // Scale by mean (mirrors R reference script exactly)
    const n = nullAD.length;
    let sum = 0;
    for (let i = 0; i < n; i++) sum += nullAD[i];
    const sFac = sum / n;
    if (sFac < 1e-12) return empP;

    const scaled = new Array(n);
    for (let i = 0; i < n; i++) scaled[i] = nullAD[i] / sFac;
    const sObs = obsAD / sFac;

    const fit = fitGenGamma(scaled);
    if (!fit.ok) return empP;

    const p = ggSurvival(sObs, fit.mu, fit.sigma, fit.Q);
    if (!isFinite(p) || p < 0 || p > 1) return empP;
    return Math.max(p, 1e-16);
  } catch {
    return empP;
  }
}

// ── Gamma distribution for KS null ───────────────────────────
//
// We fit Gamma(shape, rate) to |null_KS| via MLE.
// Parameters: [log(shape), log(rate)] — log-transform for positivity.
// Then p-value = P(Gamma > |obs_KS|) = upper tail.

function gammaNegLL(logPar, data) {
  const shape = Math.exp(logPar[0]);
  const rate  = Math.exp(logPar[1]);
  if (!isFinite(shape) || !isFinite(rate)) return 1e20;
  const js = _jstat();
  let ll = 0;
  for (let i = 0; i < data.length; i++) {
    const v = js.gamma.pdf(data[i], shape, 1 / rate);
    if (v <= 0 || !isFinite(v)) return 1e20;
    ll += Math.log(v);
  }
  return -ll;
}

/**
 * Fit Gamma(shape, rate) to positive data.
 * @param {number[]} data  positive values
 * @returns {{ shape, rate, ok }}
 */
export function fitGamma(data) {
  const pos = data.filter(v => v > 0);
  if (pos.length < 5) return { ok: false };

  // Method-of-moments seed
  let m1 = 0, m2 = 0;
  for (const v of pos) { m1 += v; m2 += v * v; }
  m1 /= pos.length; m2 /= pos.length;
  const vr = m2 - m1 * m1;
  if (vr <= 0) return { ok: false };
  const shape0 = (m1 * m1) / vr;
  const rate0  = m1 / vr;

  const result = nelderMead(
    par => gammaNegLL(par, pos),
    [Math.log(shape0), Math.log(rate0)],
    { maxIter: 1000, tol: 1e-10 }
  );

  const shape = Math.exp(result.x[0]);
  const rate  = Math.exp(result.x[1]);

  if (!isFinite(shape) || !isFinite(rate) || shape <= 0 || rate <= 0)
    return { ok: false };

  return { shape, rate, ok: true };
}

/**
 * Compute two-sided p-value for KS from fitted Gamma.
 * P(|X| >= |obs|) using upper tail of Gamma.
 *
 * @param {number}       obsKS    observed KS statistic
 * @param {Float64Array} nullKS   null distribution
 * @param {number}       empP     fallback
 * @returns {number}
 */
export function pvalGamma(obsKS, nullKS, empP) {
  try {
    const js = _jstat();
    const n = nullKS.length;
    const absNull = new Array(n);
    for (let i = 0; i < n; i++) absNull[i] = Math.abs(nullKS[i]);

    const fit = fitGamma(absNull);
    if (!fit.ok) return empP;

    // Upper tail: P(Gamma > |obs|)
    // jStat.gamma.cdf(x, shape, scale) where scale = 1/rate
    const p = 1 - js.gamma.cdf(Math.abs(obsKS), fit.shape, 1 / fit.rate);
    if (!isFinite(p) || p < 0 || p > 1) return empP;
    return Math.max(p, 1e-16);
  } catch {
    return empP;
  }
}
