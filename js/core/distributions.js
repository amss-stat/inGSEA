// ═══════════════════════════════════════════════════════════
//  core/distributions.js
//  Parametric null-distribution fitting:
//    • KS statistic  → Gamma (shape k, rate 1/θ)
//    • AD statistic  → Generalised Gamma (μ, σ, Q)  [flexsurv param.]
//
//  All fitting is done in pure JS (no WASM, no R) using
//  method-of-moments + simple MLE Newton refinement.
// ═══════════════════════════════════════════════════════════
'use strict';

// ────────────────────────────────────────────────────────────
//  Lower incomplete gamma function (regularised) P(a,x)
//  via series expansion (good for x < a+1) and
//  continued fraction (good for x >= a+1).
//  Lanczos approximation for log-Γ.
// ────────────────────────────────────────────────────────────

const LANCZOS_G = 7;
const LANCZOS_C = [
  0.99999999999980993,
  676.5203681218851,
  -1259.1392167224028,
  771.32342877765313,
  -176.61502916214059,
  12.507343278686905,
  -0.13857109526572012,
  9.9843695780195716e-6,
  1.5056327351493116e-7
];

export function logGamma(z) {
  if (z < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - logGamma(1 - z);
  z -= 1;
  let x = LANCZOS_C[0];
  for (let i = 1; i < LANCZOS_G + 2; i++) x += LANCZOS_C[i] / (z + i);
  const t = z + LANCZOS_G + 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}

/** Regularised lower incomplete gamma P(a, x) */
export function gammaIncP(a, x) {
  if (x < 0 || a <= 0) return 0;
  if (x === 0) return 0;
  if (x < a + 1) return _gammaSeries(a, x);
  return 1 - _gammaCF(a, x);
}

/** Regularised upper incomplete gamma Q(a, x) = 1 − P(a, x) */
export function gammaIncQ(a, x) { return 1 - gammaIncP(a, x); }

function _gammaSeries(a, x) {
  const lnGamA = logGamma(a);
  let ap = a, sum = 1 / a, del = 1 / a;
  for (let n = 0; n < 200; n++) {
    ap++;
    del *= x / ap;
    sum += del;
    if (Math.abs(del) < Math.abs(sum) * 3e-8) break;
  }
  return sum * Math.exp(-x + a * Math.log(x) - lnGamA);
}

function _gammaCF(a, x) {
  const lnGamA = logGamma(a);
  let b = x + 1 - a, c = 1e30, d = 1 / b, h = d;
  for (let n = 1; n <= 200; n++) {
    const an = -n * (n - a);
    b += 2;
    d = an * d + b; if (Math.abs(d) < 1e-30) d = 1e-30;
    c = b + an / c; if (Math.abs(c) < 1e-30) c = 1e-30;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 3e-8) break;
  }
  return Math.exp(-x + a * Math.log(x) - lnGamA) * h;
}

// ────────────────────────────────────────────────────────────
//  Gamma distribution fitting (for |KS| null)
//  We fit to |null_ks| values using method-of-moments then
//  one round of MLE Newton-Raphson (digamma approximation).
// ────────────────────────────────────────────────────────────

/** digamma approximation (Bernardo 1976) */
function digamma(x) {
  if (x <= 0) return NaN;
  if (x < 6) return digamma(x + 1) - 1 / x;
  return Math.log(x) - 1 / (2 * x) - 1 / (12 * x * x) + 1 / (120 * x ** 4);
}

/** trigamma approximation */
function trigamma(x) {
  if (x <= 0) return NaN;
  if (x < 6) return trigamma(x + 1) + 1 / (x * x);
  return 1 / x + 1 / (2 * x * x) + 1 / (6 * x ** 3) - 1 / (30 * x ** 5);
}

/**
 * Fit Gamma(shape, rate) to positive data.
 * Returns { shape, rate, ok }.
 */
export function fitGamma(data) {
  const n = data.length;
  if (n < 3) return { shape: 1, rate: 1, ok: false };
  let m1 = 0, m2 = 0;
  for (const v of data) { m1 += v; m2 += v * v; }
  m1 /= n; m2 /= n;
  const vr = m2 - m1 * m1;
  if (vr <= 0) return { shape: 1, rate: m1, ok: false };

  // MOM initial estimate
  let k = (m1 * m1) / vr;
  const theta = m1 / k;  // scale

  // 1 round of Newton for log-likelihood MLE
  const lnBar = data.reduce((s, v) => s + Math.log(v + 1e-15), 0) / n;
  const A = Math.log(m1) - lnBar;
  // Newton: k_new = k_old - f/f'
  for (let iter = 0; iter < 20; iter++) {
    const dk = (Math.log(k) - digamma(k) - A) / (1 / k - trigamma(k));
    k -= dk;
    if (k <= 0) { k = 1e-3; break; }
    if (Math.abs(dk) < 1e-7) break;
  }
  const rate = k / m1;  // rate = 1/scale
  return { shape: Math.max(k, 1e-4), rate: Math.max(rate, 1e-9), ok: true };
}

/**
 * Gamma CDF upper tail: P(X > x | shape, rate)
 */
export function gammaSurvival(x, shape, rate) {
  if (x <= 0) return 1;
  return gammaIncQ(shape, x * rate);
}

// ────────────────────────────────────────────────────────────
//  Generalised Gamma (flexsurv parameterisation)
//
//  GG(μ, σ, Q):
//    If Q ≠ 0:  let w = (log(t) − μ)/σ
//               u = |Q|^{-2} · exp(Q·w·sign(Q))
//               X ~ Gamma(|Q|^{-2}, 1)
//    If Q = 0:  log-normal(μ, σ)
//    If Q < 0:  reverse direction (upper tail)
//
//  We fit via:
//    1. Scale data by mean (s_fac) — matches R code exactly
//    2. log-transform → fit Normal(μ, σ²) moments (log-normal seed)
//    3. Gradient-descent MLE refining (μ, σ, Q)
// ────────────────────────────────────────────────────────────

/**
 * GG log-pdf given params.
 * Mirrors flexsurv::dgengamma with log=TRUE.
 */
function ggLogPdf(t, mu, sigma, Q) {
  if (t <= 0) return -Infinity;
  const lnt = Math.log(t);
  if (Math.abs(Q) < 1e-8) {
    // log-normal limit
    const w = (lnt - mu) / sigma;
    return -Math.log(sigma * t) - 0.5 * w * w - 0.5 * Math.log(2 * Math.PI);
  }
  const w = (lnt - mu) / sigma;
  const q2inv = 1 / (Q * Q);
  const u = q2inv * Math.exp(Q * w);
  return -Math.log(sigma * t) + logGamma(q2inv) * 0 // offset
    + (q2inv * Math.log(q2inv) - logGamma(q2inv))
    + q2inv * Q * w - u
    - Math.log(Math.abs(Q));
}

/**
 * Full log-pdf (correctly normalised).
 */
function ggLogPdfFull(t, mu, sigma, Q) {
  if (t <= 0) return -Infinity;
  if (Math.abs(Q) < 1e-8) {
    // log-normal
    const w = (Math.log(t) - mu) / sigma;
    return -Math.log(t) - Math.log(sigma) - 0.5 * Math.log(2 * Math.PI) - 0.5 * w * w;
  }
  const w = (Math.log(t) - mu) / sigma;
  const q2inv = 1 / (Q * Q);
  const u = q2inv * Math.exp(Q * w);
  // log-pdf = -log(sigma*t) + (q2inv)*Q*w - u + q2inv*log(q2inv) - logGamma(q2inv) - log|Q|
  return -Math.log(sigma) - Math.log(t)
    + q2inv * Q * w
    - u
    + q2inv * Math.log(q2inv)
    - logGamma(q2inv)
    - Math.log(Math.abs(Q));
}

/**
 * GG survival function P(T > x) = P(T > x | μ, σ, Q).
 */
export function ggSurvival(x, mu, sigma, Q) {
  if (x <= 0) return 1;
  if (Math.abs(Q) < 1e-8) {
    // log-normal: Φ(-z)
    const z = (Math.log(x) - mu) / sigma;
    return 0.5 * erfc(z / Math.SQRT2);
  }
  const w = (Math.log(x) - mu) / sigma;
  const q2inv = 1 / (Q * Q);
  const u = q2inv * Math.exp(Q * w);
  if (Q > 0) {
    // upper tail = Q(q2inv, u)  (upper regularised gamma)
    return gammaIncQ(q2inv, u);
  } else {
    // Q < 0: flip
    return gammaIncP(q2inv, u);
  }
}

/** Complementary error function */
function erfc(x) {
  // Abramowitz & Stegun 7.1.26
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  const y = 1 - poly * Math.exp(-x * x);
  return x >= 0 ? 1 - y : 1 + y;
}

/**
 * Fit Generalised Gamma to positive data.
 *
 * Strategy (mirrors R's flexsurvreg(Surv(s_data) ~ 1, dist='gengamma')):
 *   1. Scale data by mean (as in the R script).
 *   2. Seed: μ = mean(log t), σ = sd(log t), Q = 0 (log-normal seed).
 *   3. L-BFGS-style gradient descent with line search.
 *
 * @returns {{ mu, sigma, Q, ok }}
 */
export function fitGenGamma(rawData) {
  const n = rawData.length;
  if (n < 5) return { mu: 0, sigma: 1, Q: 0, ok: false };

  // Scale by mean
  const sFac = rawData.reduce((s, v) => s + v, 0) / n;
  const data = rawData.map(v => v / (sFac + 1e-15));

  // Log-transform stats for seed
  const logD = data.map(v => Math.log(Math.max(v, 1e-15)));
  const muSeed = logD.reduce((s, v) => s + v, 0) / n;
  const varSeed = logD.reduce((s, v) => s + (v - muSeed) ** 2, 0) / (n - 1);
  const sigSeed = Math.sqrt(Math.max(varSeed, 1e-6));

  // Negative log-likelihood
  const negLogLik = (mu, sigma, Q) => {
    if (sigma <= 0) return Infinity;
    let ll = 0;
    for (const t of data) ll += ggLogPdfFull(t, mu, sigma, Q);
    return -ll;
  };

  // Numerical gradient
  const grad = (mu, sigma, Q) => {
    const eps = 1e-5;
    const f0 = negLogLik(mu, sigma, Q);
    return [
      (negLogLik(mu + eps, sigma, Q) - f0) / eps,
      (negLogLik(mu, sigma + eps, Q) - f0) / eps,
      (negLogLik(mu, sigma, Q + eps) - f0) / eps
    ];
  };

  let mu = muSeed, sigma = sigSeed, Q = 0.1;
  let f = negLogLik(mu, sigma, Q);

  // Simple gradient descent with adaptive step
  let lr = 0.02;
  for (let iter = 0; iter < 400; iter++) {
    const [gMu, gSig, gQ] = grad(mu, sigma, Q);
    const norm = Math.sqrt(gMu ** 2 + gSig ** 2 + gQ ** 2) + 1e-12;

    const muNew    = mu    - lr * gMu  / norm;
    const sigmaNew = Math.max(sigma - lr * gSig / norm, 1e-4);
    const QNew     = Q     - lr * gQ   / norm;

    const fNew = negLogLik(muNew, sigmaNew, QNew);
    if (fNew < f) {
      mu = muNew; sigma = sigmaNew; Q = QNew; f = fNew;
      lr = Math.min(lr * 1.05, 0.5);
    } else {
      lr *= 0.5;
      if (lr < 1e-9) break;
    }
  }

  const ok = isFinite(f) && sigma > 1e-5;
  return { mu, sigma, Q, ok, sFac };
}

/**
 * Compute p-value from fitted GG model.
 * Scales observation by same sFac.
 * Returns empirical p if fitting failed.
 */
export function pvalFromGG(obsAD, ggFit, empiricalP) {
  if (!ggFit.ok) return empiricalP;
  const sObs = obsAD / (ggFit.sFac + 1e-15);
  try {
    const p = ggSurvival(sObs, ggFit.mu, ggFit.sigma, ggFit.Q);
    if (!isFinite(p) || p < 0 || p > 1) return empiricalP;
    return Math.max(p, 1e-16);
  } catch {
    return empiricalP;
  }
}

/**
 * Compute p-value from fitted Gamma model.
 * Uses two-sided: P(|X| >= |obs|).
 */
export function pvalFromGamma(obsKS, gammaFit, empiricalP) {
  if (!gammaFit.ok) return empiricalP;
  const p = gammaSurvival(Math.abs(obsKS), gammaFit.shape, gammaFit.rate);
  if (!isFinite(p) || p < 0 || p > 1) return empiricalP;
  return Math.max(p, 1e-16);
}
