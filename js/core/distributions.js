// ═══════════════════════════════════════════════════════════
//  core/distributions.js  ·  v2.4
//
//  Uses jStat for:
//    jStat.gammainc(x, a)  — regularised lower incomplete gamma P(a,x)
//    jStat.lngamma(a)      — log Γ(a)
//    jStat.normal.cdf(z,0,1) — standard normal CDF
//
//  Changes from previous version:
//  • gammaNegLL: replaced jStat.gamma.pdf (which can underflow)
//    with direct log-likelihood using jStat.lngamma — numerically stable
//  • ggLogPdf: verified against flexsurv source, added overflow guard on u
//  • nelderMead: added function call counter guard for robustness
// ═══════════════════════════════════════════════════════════
'use strict';

function _js() {
  if (typeof jStat === 'undefined')
    throw new Error('jStat not loaded');
  return jStat;
}

// ── Nelder-Mead (unchanged) ───────────────────────────────────
export function nelderMead(f, x0, opts = {}) {
  const {
    maxIter  = 3000,
    maxCalls = 50000,
    tol      = 1e-10,
    alpha    = 1.0,
    gamma    = 2.0,
    rho      = 0.5,
    sigma    = 0.5
  } = opts;

  const n = x0.length;
  let calls = 0;
  const _f = x => { calls++; return f(x); };

  let s  = [x0.slice()];
  for (let i = 0; i < n; i++) {
    const v = x0.slice();
    v[i] += Math.abs(v[i]) > 1e-8 ? 0.05 * Math.abs(v[i]) : 0.00025;
    s.push(v);
  }
  let fv = s.map(_f);

  for (let iter = 0; iter < maxIter && calls < maxCalls; iter++) {
    const idx = Array.from({ length: n + 1 }, (_, i) => i);
    idx.sort((a, b) => fv[a] - fv[b]);
    s  = idx.map(i => s[i]);
    fv = idx.map(i => fv[i]);

    if (fv[n] - fv[0] < tol) break;

    const c = new Array(n).fill(0);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++) c[j] += s[i][j] / n;

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

// ── GG log-pdf ────────────────────────────────────────────────
// log f(t; μ, σ, Q):
//   w = (log t − μ) / σ
//   k = 1/Q²
//   u = k · exp(Q·w)
//   log f = log|Q| + k·log(k) − logΓ(k) − log σ − log t + k·Q·w − u
//
// Q → 0 limit: log-Normal(μ, σ)

function ggLogPdf(t, mu, sigma, Q) {
  if (t <= 0 || sigma <= 0) return -Infinity;
  const js  = _js();
  const lnt = Math.log(t);

  if (Math.abs(Q) < 1e-8) {
    const z = (lnt - mu) / sigma;
    return -lnt - Math.log(sigma)
           - 0.5 * Math.log(2 * Math.PI)
           - 0.5 * z * z;
  }

  const k  = 1.0 / (Q * Q);
  const w  = (lnt - mu) / sigma;
  const Qw = Q * w;
  if (Qw > 700) return -Infinity;   // exp(Qw) would overflow → density ≈ 0
  const u = k * Math.exp(Qw);
  if (!isFinite(u) || u < 0) return -Infinity;

  return Math.log(Math.abs(Q))
       + k * Math.log(k)
       - js.gammaln(k)          // ← was js.lngamma(k)  — FIXED
       - Math.log(sigma)
       - lnt
       + k * Qw
       - u;
}

function ggNegLL(par, data) {
  const mu    = par[0];
  const sigma = Math.exp(par[1]);
  const Q     = par[2];
  if (sigma < 1e-8 || sigma > 100) return 1e15;
  let ll = 0;
  for (let i = 0; i < data.length; i++) {
    const v = ggLogPdf(data[i], mu, sigma, Q);
    if (!isFinite(v)) return 1e15;
    ll += v;
  }
  return isFinite(ll) ? -ll : 1e15;
}

// ── GG survival: P(T > x | μ, σ, Q) ──────────────────────────
// Uses jStat.gammainc(u, k) = P(k, u) = lower regularised gamma
//   Q > 0: S(x) = 1 − P(k, u)
//   Q < 0: S(x) = P(k, u)
//   Q = 0: log-normal upper tail

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

  const lowerP = js.gammainc(u, k);
  const p      = Q > 0 ? 1 - lowerP : lowerP;
  return Math.max(0, Math.min(1, p));
}

// ── Fit GG via Nelder-Mead MLE ────────────────────────────────
export function fitGenGamma(data) {
  const pos = data.filter(v => v > 0);
  if (pos.length < 10) return { ok: false };

  const logD = pos.map(v => Math.log(v));
  const n    = pos.length;
  const muS  = logD.reduce((s, v) => s + v, 0) / n;
  const varS = logD.reduce((s, v) => s + (v - muS) ** 2, 0) / (n - 1);
  const sigS = Math.max(Math.sqrt(varS), 1e-4);

  // Three starts to avoid local minima (Q < 0, Q ≈ 0, Q > 0)
  const starts = [
    [muS, Math.log(sigS), -1.0],
    [muS, Math.log(sigS),  0.1],
    [muS, Math.log(sigS),  1.0],
  ];

  let best = null;
  for (const start of starts) {
    const res = nelderMead(
      par => ggNegLL(par, pos),
      start,
      { maxIter: 3000, tol: 1e-10 }
    );
    if (!best || res.fval < best.fval) best = res;
  }

  if (!best || !isFinite(best.fval)) return { ok: false };

  const mu    = best.x[0];
  const sigma = Math.exp(best.x[1]);
  const Q     = best.x[2];

  if (sigma < 1e-8 || sigma > 50) return { ok: false };
  return { mu, sigma, Q, ok: true };
}

// ── p-value for AD from GG fit ────────────────────────────────
export function pvalGG(obsAD, nullAD, empP) {
  try {
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
    return isFinite(p) && p >= 0 && p <= 1
      ? Math.max(p, 1e-16) : empP;
  } catch { return empP; }
}

// ── Gamma log-likelihood (numerically stable) ─────────────────
// log f(x; k, r) = (k−1)·log x + k·log r − r·x − logΓ(k)
// Parameters: [log k, log r]  (log-transform ensures positivity)

function gammaNegLL(logPar, data) {
  const k = Math.exp(logPar[0]);
  const r = Math.exp(logPar[1]);
  if (!isFinite(k) || !isFinite(r) || k < 1e-6 || r < 1e-6) return 1e15;
  const js     = _js();
  const lnGamK = js.gammaln(k);   // ← was js.lngamma(k)  — FIXED
  const lnR    = logPar[1];
  let ll = 0;
  for (let i = 0; i < data.length; i++) {
    const x = data[i];
    if (x <= 0) return 1e15;
    ll += (k - 1) * Math.log(x) + k * lnR - r * x - lnGamK;
  }
  return isFinite(ll) ? -ll : 1e15;
}

// ── Fit Gamma to positive data ────────────────────────────────
export function fitGamma(data) {
  const pos = data.filter(v => v > 0);
  if (pos.length < 10) return { ok: false };

  let m1 = 0, m2 = 0;
  for (const v of pos) { m1 += v; m2 += v * v; }
  m1 /= pos.length; m2 /= pos.length;
  const vr = m2 - m1 * m1;
  if (vr < 1e-14) return { ok: false };

  const k0 = (m1 * m1) / vr;
  const r0 = m1 / vr;

  const res = nelderMead(
    par => gammaNegLL(par, pos),
    [Math.log(k0), Math.log(r0)],
    { maxIter: 1000, tol: 1e-10 }
  );

  const k = Math.exp(res.x[0]);
  const r = Math.exp(res.x[1]);
  if (!isFinite(k) || !isFinite(r) || k <= 0 || r <= 0)
    return { ok: false };
  return { shape: k, rate: r, ok: true };
}

// ── p-value for KS from Gamma fit ────────────────────────────
// Two-sided: P(|X| ≥ |obs|) = upper tail of Gamma fitted to |null_KS|
// jStat.gamma.cdf(x, shape, scale) where scale = 1/rate

export function pvalGamma(obsKS, nullKS, empP) {
  try {
    const js = _js();
    const n  = nullKS.length;
    const absNull = new Array(n);
    for (let i = 0; i < n; i++) absNull[i] = Math.abs(nullKS[i]);

    const fit = fitGamma(absNull);
    if (!fit.ok) return empP;

    const p = 1 - js.gamma.cdf(Math.abs(obsKS), fit.shape, 1 / fit.rate);
    return isFinite(p) && p >= 0 && p <= 1
      ? Math.max(p, 1e-16) : empP;
  } catch { return empP; }
}
