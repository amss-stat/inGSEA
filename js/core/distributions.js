// ═══════════════════════════════════════════════════════════
//  core/distributions.js  ·  v2.9 (Performance Optimized)
//
//  Key Optimizations:
//  • Pre-calculates Math.log(data) to avoid O(Iter * N) redundant logs.
//  • Subsamples Null Distribution to 400 points (sufficient for GG fit).
//  • Vectorized-style NegLL: constant terms moved out of loops.
//  • Reduced starting points and max iterations for Nelder-Mead.
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
    maxIter  = 1000, // 降低默认迭代上限
    maxCalls = 5000,
    tol      = 1e-8, // 调整容差到合理范围
    alpha    = 1.0,
    gamma    = 2.0,
    rho      = 0.5,
    sigma    = 0.5
  } = opts;

  const n = x0.length;
  let calls = 0;
  const _f = x => { calls++; return f(x); };

  let s = [x0.slice()];
  for (let i = 0; i < n; i++) {
    const v = x0.slice();
    v[i] += (Math.abs(v[i]) > 1e-8) ? 0.05 * Math.abs(v[i]) : 0.00025;
    s.push(v);
  }
  let fv = s.map(_f);

  for (let iter = 0; iter < maxIter && calls < maxCalls; iter++) {
    const idx = Array.from({ length: n + 1 }, (_, i) => i);
    idx.sort((a, b) => fv[a] - fv[b]);
    s = idx.map(i => s[i]);
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

// ── Optimized Generalized Gamma Likelihood ───────────────────

/**
 * 优化后的负对数似然函数
 * @param {Array} par - [mu, logSigma, Q]
 * @param {Array} logData - 预计算好的对数数组
 */
function ggNegLL_Optimized(par, logData) {
  const mu    = par[0];
  const sigma = Math.exp(par[1]);
  const Q     = par[2];
  const n     = logData.length;
  const js    = _js();

  if (sigma < 1e-7 || sigma > 200) return 1e15;

  // 情况 A: 接近对数正态 (Q -> 0)
  if (Math.abs(Q) < 1e-6) {
    let sumZ2 = 0;
    for (let i = 0; i < n; i++) {
      const z = (logData[i] - mu) / sigma;
      sumZ2 += z * z;
    }
    return n * (Math.log(sigma) + 0.918938533) + 0.5 * sumZ2;
  }

  // 情况 B: 广义伽马
  const k = 1.0 / (Q * Q);
  const invSigma = 1.0 / sigma;
  const kQ = k * Q;
  
  // 提取循环外的常数项计算
  // logF = log|Q| + k*log(k) - logGamma(k) - log(sigma) - log(t) + k*Q*w - u
  const constTerms = Math.log(Math.abs(Q)) + k * Math.log(k) - js.gammaln(k) - Math.log(sigma);
  
  let sumDynamic = 0;
  for (let i = 0; i < n; i++) {
    const lnt = logData[i];
    const w   = (lnt - mu) * invSigma;
    const Qw  = Q * w;
    if (Qw > 700) return 1e15; // 溢出保护
    const u   = k * Math.exp(Qw);
    sumDynamic += lnt - kQ * w + u;
  }

  return -(n * constTerms - sumDynamic);
}

// ── Survival Function ────────────────────────────────────────

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

  const lowerP = js.lowRegGamma(k, u);
  const p = Q > 0 ? 1 - lowerP : lowerP;
  return Math.max(0, Math.min(1, p));
}

// ── Fitting Logic ────────────────────────────────────────────

export function fitGenGamma(data) {
  // 1. 过滤并下采样
  const pos = data.filter(v => v > 1e-12);
  if (pos.length < 10) return { ok: false };

  let fitData = pos;
  if (pos.length > 1000) {
    fitData = [];
    const step = pos.length / 1000;
    for (let i = 0; i < 1000; i++) {
      fitData.push(pos[Math.floor(i * step)]);
    }
  }

  // 2. 关键优化：预计算对数，避免在优化器迭代中重复计算
  const logData = fitData.map(v => Math.log(v));
  const n = logData.length;

  const muS  = logData.reduce((s, v) => s + v, 0) / n;
  const varS = logData.reduce((s, v) => s + (v - muS) ** 2, 0) / (n - 1);
  const sigS = Math.max(Math.sqrt(varS), 1e-4);

  // 3. 减少起点数量，覆盖正偏和负偏即可
  const starts = [
    [muS, Math.log(sigS), -1.0],
    [muS, Math.log(sigS),  0.8]
  ];

  let best = null;
  for (const start of starts) {
    const res = nelderMead(par => ggNegLL_Optimized(par, logData), start, {
      maxIter: 1000,
      tol: 1e-8
    });
    if (best === null || res.fval < best.fval) best = res;
  }

  if (!best || !isFinite(best.fval)) return { ok: false };

  const mu    = best.x[0];
  const sigma = Math.exp(best.x[1]);
  const Q     = best.x[2];

  if (sigma < 1e-7 || sigma > 200) return { ok: false };
  return { mu, sigma, Q, ok: true };
}

/**
 * 接口函数：计算参数化 p 值
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
  } catch (err) {
    console.warn('GG Fit Error:', err.message);
    return { p: empP, fitted: false };
  }
}
