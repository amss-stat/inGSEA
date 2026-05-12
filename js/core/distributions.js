// ═══════════════════════════════════════════════════════════
//  core/distributions.js  ·  v2.7
//
//  Uses jStat for:
//    jStat.gammainc(x, a)  — regularised lower incomplete gamma P(a,x)
//    jStat.lngamma(a)      — log Γ(a)
//    jStat.normal.cdf(z,0,1) — standard normal CDF
//
//  Changes from v2.6:
//  • Removed Gamma distribution fitting (pvalGamma / fitGamma)
//    KS now uses permutation p-values exclusively
//  • pvalGG returns { p, fitted } structured result
//    instead of sentinel-value fallback
//  • fitGenGamma: relaxed sigma guard (50 → 200),
//    added more Nelder-Mead starts, increased iterations
//  • Positive-data filter uses 1e-12 threshold instead of 0
// ═══════════════════════════════════════════════════════════
'use strict';

function _js() {
  if (typeof jStat === 'undefined')
    throw new Error('jStat not loaded');
  return jStat;
}

// ── Nelder-Mead simplex optimizer ────────────────────────────
export function nelderMead(f, x0, opts = {}) {
  const {
    maxIter  = 5000,
    maxCalls = 80000,
    tol      = 1e-12,
    alpha    = 1.0,
    gamma    = 2.0,
    rho      = 0.5,
    sigma    = 0.5
  } = opts;

  const n = x0.length;
  let calls = 0;
  const _f = x => { calls++; return f(x); };

  // Initial simplex
  let s = [x0.slice()];
  for (let i = 0; i < n; i++) {
    const v = x0.slice();
    v[i] += (Math.abs(v[i]) > 1e-8) ? 0.05 * Math.abs(v[i]) : 0.00025;
    s.push(v);
  }
  let fv = s.map(_f);

  for (let iter = 0; iter < maxIter && calls < maxCalls; iter++) {
    // Sort
    const idx = Array.from({ length: n + 1 }, (_, i) => i);
    idx.sort((a, b) => fv[a] - fv[b]);
    s = idx.map(i => s[i]);
    fv = idx.map(i => fv[i]);

    if (fv[n] - fv[0] < tol) break;

    // Centroid of all but worst
    const c = new Array(n).fill(0);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++) c[j] += s[i][j] / n;

    // Reflect
    const xr = c.map((ci, j) => ci + alpha * (ci - s[n][j]));
    const fr = _f(xr);

    if (fr < fv[0]) {
      const xe = c.map((ci, j) => ci + gamma * (xr[j] - ci));
      const fe = _f(xe);
      if (fe < fr) { s[n] = xe; fv[n] = fe; }
      else          { s[n] = xr; fv[n] = fr; }
    } else if (fr < fv[n - 1]) {
      s[n] = xr; fv[n] = fr;
    } else {
      const xc = c.map((ci, j) => ci + rho * (s[n][j] - ci));
      const fc = _f(xc);
      if (fc < fv[n]) {
        s[n] = xc; fv[n] = fc;
      } else {
        for (let i = 1; i <= n; i++) {
          s[i] = s[i].map((v, j) => s[0][j] + sigma * (v - s[0][j]));
          fv[i] = _f(s[i]);
        }
      }
    }
  }

  return { x: s[0], fval: fv[0] };
}

// ── Generalised Gamma (flexsurv parameterisation) ─────────────
//
//  log f(t; μ, σ, Q):
//    w = (log t − μ) / σ
//    k = 1/Q²
//    u = k · exp(Q·w)
//    log f = log|Q| + k·log k − logΓ(k) − log σ − log t + k·Q·w − u
//
//  Q → 0: log-normal(μ, σ)

function ggLogPdf(t, mu, sigma, Q) {
  if (t <= 0 || sigma <= 0) return -Infinity;
  const lnt = Math.log(t);
  const js  = _js();

  if (Math.abs(Q) < 1e-8) {
    const z = (lnt - mu) / sigma;
    return -lnt - Math.log(sigma) - 0.5 * Math.log(2 * Math.PI) - 0.5 * z * z;
  }

  const k = 1.0 / (Q * Q);
  const w = (lnt - mu) / sigma;
  // Guard against overflow: Q*w can be large
  const Qw = Q * w;
  if (Qw > 700) return -Infinity;  // exp(Qw) would overflow
  const u = k * Math.exp(Qw);
  if (!isFinite(u) || u < 0) return -Infinity;

  return Math.log(Math.abs(Q))
    + k * Math.log(k)
    - js.lngamma(k)
    - Math.log(sigma)
    - lnt
    + k * Qw
    - u;
}

function ggNegLL(par, data) {
  const mu    = par[0];
  const sigma = Math.exp(par[1]);   // log-transform: sigma > 0 always
  const Q     = par[2];
  if (sigma < 1e-8 || sigma > 200) return 1e15;
  let ll = 0;
  for (let i = 0; i < data.length; i++) {
    const v = ggLogPdf(data[i], mu, sigma, Q);
    if (!isFinite(v)) return 1e15;
    ll += v;
  }
  return isFinite(ll) ? -ll : 1e15;
}

/**
 * GG survival: P(T > x | μ, σ, Q).
 * Uses jStat.gammainc for the regularised incomplete gamma.
 *
 * jStat.gammainc(x, a) = P(a,x) = lower regularised gamma
 *   Q > 0: P(T > x) = 1 − P(k, u)  [upper tail of Gamma(k,1)]
 *   Q < 0: P(T > x) = P(k, u)       [lower tail, direction flips]
 *   Q = 0: log-normal upper tail
 */
export function ggSurvival(x, mu, sigma, Q) {
  if (x <= 0) return 1;
  const js = _js();

  if (Math.abs(Q) < 1e-8) {
    const z = (Math.log(x) - mu) / sigma;
    return 1 - js.normal.cdf(z, 0, 1);
  }

  const k  = 1.0 / (Q * Q);
  const w  = (Math.log(x) - mu) / sigma;
  const Qw = Q * w;
  if (Qw > 700) return Q > 0 ? 0 : 1;
  const u = k * Math.exp(Qw);
  if (!isFinite(u) || u < 0) return Q > 0 ? 1 : 0;

  // jStat.gammainc(u, k) = P(k, u) = lower regularised gamma
  const lowerP = js.gammainc(u, k);
  const p = Q > 0 ? 1 - lowerP : lowerP;
  return Math.max(0, Math.min(1, p));
}

/**
 * Fit GG to positive data using Nelder-Mead MLE.
 * Returns { mu, sigma, Q, ok }.
 */
export function fitGenGamma(data) {
  const pos = data.filter(v => v > 1e-12);
  if (pos.length < 10) return { ok: false };

  const logD = pos.map(v => Math.log(v));
  const n    = pos.length;
  const muS  = logD.reduce((s, v) => s + v, 0) / n;
  const varS = logD.reduce((s, v) => s + (v - muS) ** 2, 0) / (n - 1);
  const sigS = Math.max(Math.sqrt(varS), 1e-4);

  // Multiple starts to avoid local minima
  const starts = [
    [muS, Math.log(sigS), -1.5],
    [muS, Math.log(sigS), -0.5],
    [muS, Math.log(sigS),  0.1],
    [muS, Math.log(sigS),  0.5],
    [muS, Math.log(sigS),  1.5],
  ];

  let best = null;
  for (const start of starts) {
    const res = nelderMead(par => ggNegLL(par, pos), start,
      { maxIter: 5000, tol: 1e-12 });
    if (best === null || res.fval < best.fval) best = res;
  }

  if (!best || !isFinite(best.fval)) return { ok: false };

  const mu    = best.x[0];
  const sigma = Math.exp(best.x[1]);
  const Q     = best.x[2];

  if (sigma < 1e-8 || sigma > 200) return { ok: false };

  return { mu, sigma, Q, ok: true };
}

/**
 * Upper-tail p-value for AD from GG fit.
 * Scales data by mean (matching R reference script).
 *
 * Returns { p, fitted }:
 *   fitted=true  → parametric p-value from GG survival
 *   fitted=false → fell back to empP
 */
export function pvalGG(obsAD, nullAD, empP) {
  try {
    const n = nullAD.length;
    let sum = 0;
    for (let i = 0; i < n; i++) sum += nullAD[i];
    const sFac = sum / n;
    if (sFac < 1e-12) return { p: empP, fitted: false };

    const scaled = new Array(n);
    for (let i = 0; i < n; i++) scaled[i] = nullAD[i] / sFac;
    const sObs = obsAD / sFac;

    const fit = fitGenGamma(scaled);
    if (!fit.ok) return { p: empP, fitted: false };

    const p = ggSurvival(sObs, fit.mu, fit.sigma, fit.Q);
    if (isFinite(p) && p >= 0 && p <= 1) {
      return { p: Math.max(p, 1e-16), fitted: true };
    }
    return { p: empP, fitted: false };
  } catch {
    return { p: empP, fitted: false };
  }
}
