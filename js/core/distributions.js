'use strict';

function _js() {
  if (typeof jStat === 'undefined')
    throw new Error('jStat not loaded');
  return jStat;
}

// ── Nelder-Mead simplex optimizer ────────────────────────────
export function nelderMead(f, x0, opts = {}) {
  const {
    maxIter  = 1000,
    maxCalls = 5000,
    tol      = 1e-6,
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

function ggNegLL_Optimized(par, logData, sumLogData) {
  const mu    = par[0];
  const sigma = Math.exp(par[1]);
  const Q     = par[2];
  const n     = logData.length;
  const js    = _js();

  if (sigma < 1e-7 || sigma > 200) return 1e15;

  if (Math.abs(Q) < 1e-6) {
    let sumZ2 = 0;
    const invSigma = 1.0 / sigma;
    for (let i = 0; i < n; i++) {
      const z = (logData[i] - mu) * invSigma;
      sumZ2 += z * z;
    }
    return n * (Math.log(sigma) + 0.9189385332) + sumLogData + 0.5 * sumZ2;
  }

  const k = 1.0 / (Q * Q);
  const invSigma = 1.0 / sigma;
  const kQ = k * Q;
  
  const term1 = n * (js.gammaln(k) + Math.log(sigma) - Math.log(Math.abs(Q)) - k * Math.log(k));
  
  let sumExpPart = 0;
  let sumWPart   = 0;
  for (let i = 0; i < n; i++) {
    const w = (logData[i] - mu) * invSigma;
    const Qw = Q * w;
    if (Qw > 700) return 1e15; 
    sumExpPart += k * Math.exp(Qw);
    sumWPart   += kQ * w;
  }

  return sumLogData + term1 - sumWPart + sumExpPart;
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


export function fitGenGamma(data) {
  const pos = data.filter(v => v > 1e-12);
  if (pos.length < 10) return { ok: false };


  const logData = pos.map(v => Math.log(v));
  const sumLogData = logData.reduce((a, b) => a + b, 0);
  
  const n = logData.length;
  const muS  = sumLogData / n;
  const varS = logData.reduce((s, v) => s + (v - muS) ** 2, 0) / (n - 1);
  const sigS = Math.max(Math.sqrt(varS), 1e-4);

  const starts = [
    [muS, Math.log(sigS), -1.0],
    [muS, Math.log(sigS),  0.8]
  ];

  let best = null;
  for (const start of starts) {
    const res = nelderMead(par => ggNegLL_Optimized(par, logData, sumLogData), start, {
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
