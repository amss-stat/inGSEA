// ═══════════════════════════════════════════════════════════
//  core/webr-bridge.js  ·  WebR integration
//
//  KEY DESIGN: we do NOT use flexsurv because its dependencies
//  (deSolve, muhaz) are not available as WebR WASM binaries.
//
//  Instead we implement GG fitting in pure base R using:
//    - optim() with L-BFGS-B (base R)
//    - pgamma(), lgamma(), dgamma() (base R)
//    - Manual GG log-likelihood, CDF via pgamma (base R)
//
//  This matches the mathematical model in flexsurv exactly
//  (same parameterisation: mu, sigma, Q) but uses only
//  base R functions that are guaranteed available in WebR.
// ═══════════════════════════════════════════════════════════
'use strict';

let _webR     = null;
let _ready    = false;
let _failed   = false;
let _initProm = null;

// ── Public API ───────────────────────────────────────────────

/**
 * Begin loading WebR in background.
 * @param {(state:string, msg?:string)=>void} onStatus
 * @returns {Promise<boolean>}
 */
export function initWebR(onStatus) {
  if (_initProm) return _initProm;
  _initProm = _doInit(onStatus);
  return _initProm;
}

export const isWebRReady  = () => _ready;
export const isWebRFailed = () => _failed;

async function _doInit(onStatus) {
  try {
    onStatus('loading');

    const { WebR } = await import(
      'https://webr.r-wasm.org/v0.4.2/webr.mjs'
    );
    _webR = new WebR();
    await _webR.init();

    onStatus('installing');

    // Define all fitting functions in pure base R — no external packages
    await _webR.evalRVoid(R_SETUP_CODE);

    _ready = true;
    onStatus('ready');
    return true;

  } catch (err) {
    _failed = true;
    onStatus('error', String(err?.message || err));
    return false;
  }
}

/**
 * Fit parametric null distributions for multiple pathways.
 * Each item: { nullKS:Float64Array, nullAD:Float64Array,
 *              obsKS:number, obsAD:number,
 *              empKS:number, empAD:number }
 * Returns: Array<{ pKS:number, pAD:number, engine:string }>
 */
export async function fitMany(items) {
  if (!_ready) {
    return items.map(it => ({
      pKS: it.empKS, pAD: it.empAD, engine: 'permutation'
    }));
  }

  try {
    // Build R source that defines data and calls our fitting function
    const rItems = items.map((it, i) => `
      list(
        null_ks = c(${_vec(it.nullKS)}),
        null_ad = c(${_vec(it.nullAD)}),
        obs_ks  = ${_num(it.obsKS)},
        obs_ad  = ${_num(it.obsAD)},
        emp_ks  = ${_num(it.empKS)},
        emp_ad  = ${_num(it.empAD)}
      )`).join(',\n');

    const rCode = `
      .igsea_fit_batch(list(${rItems}))
    `;

    const result = await _webR.evalR(rCode);
    const js     = await result.toJs();

    // js is a list of lists; parse out p_ks and p_ad
    return js.values.map((entry, i) => {
      const pKS = _extractVal(entry, 'p_ks', items[i].empKS);
      const pAD = _extractVal(entry, 'p_ad', items[i].empAD);
      return {
        pKS: _guardP(pKS, items[i].empKS),
        pAD: _guardP(pAD, items[i].empAD),
        engine: 'parametric'
      };
    });

  } catch (err) {
    console.warn('[iGSEA] WebR fitMany error:', err);
    return items.map(it => ({
      pKS: it.empKS, pAD: it.empAD, engine: 'permutation-fallback'
    }));
  }
}

// ── Helpers ──────────────────────────────────────────────────

function _vec(arr) {
  const parts = new Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    parts[i] = arr[i].toPrecision(8);
  }
  return parts.join(',');
}

function _num(v) {
  if (!isFinite(v)) return 'NaN';
  return v.toPrecision(12);
}

function _guardP(p, fallback) {
  if (p == null || !isFinite(p) || p < 0 || p > 1) return fallback;
  return Math.max(p, 1e-16);
}

function _extractVal(entry, name, fallback) {
  try {
    const idx = entry.names.indexOf(name);
    if (idx < 0) return fallback;
    const v = entry.values[idx].values[0];
    return isFinite(v) ? v : fallback;
  } catch { return fallback; }
}

// ── Pure base R code for GG/Gamma fitting ────────────────────
// This is injected once at init time.
// Mirrors the flexsurv GG parameterisation exactly:
//   w = (log(t) - mu) / sigma
//   If Q != 0:  u = Q^{-2} * exp(Q*w)
//               X ~ Gamma(Q^{-2}, 1)
//   If Q == 0:  log-normal(mu, sigma)
//
const R_SETUP_CODE = `
# ─────────────────────────────────────────────────────────
# Generalised Gamma: log-pdf, survival, fitting via optim()
# Parameterisation matches flexsurv (mu, sigma, Q)
# ─────────────────────────────────────────────────────────

.gg_logpdf <- function(t, mu, sigma, Q) {
  # Returns log f(t | mu, sigma, Q)
  lnt <- log(t)
  w   <- (lnt - mu) / sigma

  if (abs(Q) < 1e-8) {
    # Log-normal limit
    return(-lnt - log(sigma) - 0.5*log(2*pi) - 0.5*w^2)
  }

  k <- 1 / (Q^2)              # shape parameter
  u <- k * exp(Q * w)         # transformed variable

  # log f = log|Q| + k*log(k) - lgamma(k) - log(sigma) - log(t)
  #         + k*Q*w - u
  log(abs(Q)) + k*log(k) - lgamma(k) - log(sigma) - lnt +
    k*Q*w - u
}

.gg_negloglik <- function(par, data) {
  mu <- par[1]; sigma <- exp(par[2]); Q <- par[3]
  if (sigma < 1e-6) return(1e20)
  ll <- sum(sapply(data, function(t) .gg_logpdf(t, mu, sigma, Q)))
  if (!is.finite(ll)) return(1e20)
  -ll
}

.gg_survival <- function(x, mu, sigma, Q) {
  # P(T > x | mu, sigma, Q)
  if (x <= 0) return(1)

  if (abs(Q) < 1e-8) {
    # Log-normal survival
    z <- (log(x) - mu) / sigma
    return(pnorm(-z))
  }

  w <- (log(x) - mu) / sigma
  k <- 1 / (Q^2)
  u <- k * exp(Q * w)

  if (Q > 0) {
    # Upper tail of Gamma(k, 1)
    return(pgamma(u, shape = k, rate = 1, lower.tail = FALSE))
  } else {
    # Q < 0: lower tail
    return(pgamma(u, shape = k, rate = 1, lower.tail = TRUE))
  }
}

.fit_gg <- function(data, emp_p) {
  # Fit GG to positive data, return upper-tail p-value at obs
  # data = scaled null distribution, obs already in data context
  tryCatch({
    ld <- log(data[data > 0])
    if (length(ld) < 5) return(emp_p)
    mu0    <- mean(ld)
    sigma0 <- sd(ld)
    if (sigma0 < 1e-8) return(emp_p)

    fit <- optim(
      par     = c(mu0, log(sigma0), 0.1),
      fn      = .gg_negloglik,
      data    = data[data > 0],
      method  = "Nelder-Mead",
      control = list(maxit = 2000, reltol = 1e-10)
    )

    mu    <- fit$par[1]
    sigma <- exp(fit$par[2])
    Q     <- fit$par[3]

    if (!is.finite(fit$value) || sigma < 1e-6) return(emp_p)
    return(c(mu = mu, sigma = sigma, Q = Q))
  }, error = function(e) {
    return(emp_p)
  })
}

# ─────────────────────────────────────────────────────────
# Gamma fitting for |KS| null
# Uses method of moments + optim MLE refinement
# ─────────────────────────────────────────────────────────

.fit_gamma_ks <- function(abs_null_ks, obs_abs_ks, emp_p) {
  tryCatch({
    x <- abs_null_ks[abs_null_ks > 0]
    if (length(x) < 5) return(emp_p)

    m  <- mean(x)
    v  <- var(x)
    if (v < 1e-12) return(emp_p)

    shape0 <- m^2 / v
    rate0  <- m / v

    # MLE refinement
    negll <- function(par) {
      shape <- exp(par[1]); rate <- exp(par[2])
      -sum(dgamma(x, shape = shape, rate = rate, log = TRUE))
    }

    fit <- optim(
      par     = c(log(shape0), log(rate0)),
      fn      = negll,
      method  = "Nelder-Mead",
      control = list(maxit = 1000, reltol = 1e-10)
    )

    shape <- exp(fit$par[1])
    rate  <- exp(fit$par[2])

    p <- pgamma(obs_abs_ks, shape = shape, rate = rate, lower.tail = FALSE)
    return(max(min(p, 1), 1e-16))
  }, error = function(e) emp_p)
}

# ─────────────────────────────────────────────────────────
# Batch fitting: single R call for all pathways
# ─────────────────────────────────────────────────────────

.igsea_fit_one <- function(item) {
  null_ks <- item$null_ks
  null_ad <- item$null_ad
  obs_ks  <- item$obs_ks
  obs_ad  <- item$obs_ad
  emp_ks  <- item$emp_ks
  emp_ad  <- item$emp_ad

  # --- KS: Gamma on |null_ks| ---
  p_ks <- .fit_gamma_ks(abs(null_ks), abs(obs_ks), emp_ks)

  # --- AD: Generalised Gamma on scaled null_ad ---
  s_fac <- mean(null_ad)
  if (s_fac < 1e-12) return(list(p_ks = emp_ks, p_ad = emp_ad))

  s_data <- null_ad / s_fac
  s_obs  <- obs_ad  / s_fac

  gg_result <- .fit_gg(s_data, emp_ad)

  if (is.numeric(gg_result) && length(gg_result) == 3) {
    mu    <- gg_result["mu"]
    sigma <- gg_result["sigma"]
    Q     <- gg_result["Q"]
    p_ad  <- .gg_survival(s_obs, mu, sigma, Q)
    p_ad  <- max(min(p_ad, 1), 1e-16)
  } else {
    p_ad <- emp_ad
  }

  list(p_ks = as.numeric(p_ks), p_ad = as.numeric(p_ad))
}

.igsea_fit_batch <- function(items) {
  lapply(items, .igsea_fit_one)
}

cat("[iGSEA] Pure base-R fitting functions loaded.\\n")
`;
