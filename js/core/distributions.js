// ═══════════════════════════════════════════════════════════
//  core/distributions.js  ·  v2.8 (Zero Dependency Edition)
//
//  Changes:
//  • Removed jStat entirely. Completely self-contained.
//  • Added built-in numerical approximations for:
//      - gammaln(x): Log-Gamma (Lanczos approximation)
//      - normalCDF(x): Standard Normal CDF (Abramowitz & Stegun)
//      - lowRegGamma(a,x): Lower Regularized Incomplete Gamma P(a,x)
// ═══════════════════════════════════════════════════════════
'use strict';

// ── 内部数学核心库 (替换 jStat) ─────────────────────────────

/** 计算 Gamma 函数的自然对数: ln(Γ(x)) (Lanczos 近似) */
function gammaln(x) {
  const cof = [
    76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
  ];
  let y = x;
  let tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) {
    y += 1;
    ser += cof[j] / y;
  }
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

/** 标准正态分布累积概率 CDF: Φ(x) */
function normalCDF(x) {
  const t = 1 / (1 + 0.2316419 * Math.abs(x));
  const d = 0.3989423 * Math.exp(-x * x / 2);
  const p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
  return x > 0 ? 1 - p : p;
}

/** 下不完全 Gamma 函数 P(a, x) = γ(a,x)/Γ(a) */
function lowRegGamma(a, x) {
  if (x < 0 || a <= 0) return 0;
  if (x === 0) return 0;

  const gln = gammaln(a);

  if (x < a + 1.0) {
    // 级数展开 (Series Representation)
    let ap = a;
    let sum = 1.0 / a;
    let del = sum;
    for (let n = 1; n <= 100; n++) {
      ap += 1.0;
      del *= x / ap;
      sum += del;
      if (Math.abs(del) < Math.abs(sum) * 1e-14) break;
    }
    return sum * Math.exp(-x + a * Math.log(x) - gln);
  } else {
    // 连分数展开 (Continued Fraction), 计算上 Gamma 然后用 1 减去
    let b = x + 1.0 - a;
    let c = 1.0 / 1e-30;
    let d = 1.0 / b;
    let h = d;
    for (let i = 1; i <= 100; i++) {
      let an = -i * (i - a);
      b += 2.0;
      d = an * d + b;
      if (Math.abs(d) < 1e-30) d = 1e-30;
      c = b + an / c;
      if (Math.abs(c) < 1e-30) c = 1e-30;
      d = 1.0 / d;
      let del = d * c;
      h *= del;
      if (Math.abs(del - 1.0) < 1e-14) break;
    }
    const upper = Math.exp(-x + a * Math.log(x) - gln) * h;
    return 1.0 - upper;
  }
}

// ── Nelder-Mead simplex optimizer ────────────────────────────
export function nelderMead(f, x0, opts = {}) {
  // ...保持原有代码不变...
  const { maxIter = 5000, maxCalls = 80000, tol = 1e-12, alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5 } = opts;
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
      else         { s[n] = xr; fv[n] = fr; }
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

// ── Generalised Gamma ─────────────────────────────
function ggLogPdf(t, mu, sigma, Q) {
  if (t <= 0 || sigma <= 0) return -Infinity;
  const lnt = Math.log(t);

  if (Math.abs(Q) < 1e-8) {
    const z = (lnt - mu) / sigma;
    return -lnt - Math.log(sigma) - 0.5 * Math.log(2 * Math.PI) - 0.5 * z * z;
  }

  const k = 1.0 / (Q * Q);
  const w = (lnt - mu) / sigma;
  const Qw = Q * w;
  if (Qw > 700) return -Infinity;  
  const u = k * Math.exp(Qw);
  if (!isFinite(u) || u < 0) return -Infinity;

  return Math.log(Math.abs(Q))
    + k * Math.log(k)
    - gammaln(k) // 使用内置方法
    - Math.log(sigma)
    - lnt
    + k * Qw
    - u;
}

function ggNegLL(par, data) {
  const mu    = par[0];
  const sigma = Math.exp(par[1]);
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

export function ggSurvival(x, mu, sigma, Q) {
  if (x <= 0) return 1;

  if (Math.abs(Q) < 1e-8) {
    const z = (Math.log(x) - mu) / sigma;
    return 1 - normalCDF(z); // 使用内置方法
  }

  const k  = 1.0 / (Q * Q);
  const w  = (Math.log(x) - mu) / sigma;
  const Qw = Q * w;
  if (Qw > 700) return Q > 0 ? 0 : 1;
  const u = k * Math.exp(Qw);
  if (!isFinite(u) || u < 0) return Q > 0 ? 1 : 0;

  const lowerP = lowRegGamma(k, u); // 使用内置方法
  const p = Q > 0 ? 1 - lowerP : lowerP;
  return Math.max(0, Math.min(1, p));
}

export function fitGenGamma(data) {
  // ...保持原有代码不变...
  const pos = data.filter(v => v > 1e-12);
  if (pos.length < 10) return { ok: false };
  const logD = pos.map(v => Math.log(v));
  const n    = pos.length;
  const muS  = logD.reduce((s, v) => s + v, 0) / n;
  const varS = logD.reduce((s, v) => s + (v - muS) ** 2, 0) / (n - 1);
  const sigS = Math.max(Math.sqrt(varS), 1e-4);
  const starts = [
    [muS, Math.log(sigS), -1.5],
    [muS, Math.log(sigS), -0.5],
    [muS, Math.log(sigS),  0.1],
    [muS, Math.log(sigS),  0.5],
    [muS, Math.log(sigS),  1.5],
  ];
  let best = null;
  for (const start of starts) {
    const res = nelderMead(par => ggNegLL(par, pos), start, { maxIter: 5000, tol: 1e-12 });
    if (best === null || res.fval < best.fval) best = res;
  }
  if (!best || !isFinite(best.fval)) return { ok: false };
  const mu    = best.x[0];
  const sigma = Math.exp(best.x[1]);
  const Q     = best.x[2];
  if (sigma < 1e-8 || sigma > 200) return { ok: false };
  return { mu, sigma, Q, ok: true };
}

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
    console.warn('Parametric fit failed:', err);
    return { p: empP, fitted: false };
  }
}
