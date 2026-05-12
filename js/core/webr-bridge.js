// ═══════════════════════════════════════════════════════════
//  core/webr-bridge.js  ·  WebR + flexsurv integration
//
//  WebR v0.4.x API reference:
//    https://docs.r-wasm.org/webr/latest/api/js/
//
//  Key design decisions:
//  1. Load WebR once at startup (background, non-blocking).
//  2. Each fitting call uses evalR() with data serialized
//     as R source strings (most reliable cross-version approach).
//  3. All R errors are caught; fall back to empirical p.
//  4. fitMany() batches all pathways in a single R call
//     to minimise JS↔R round-trip overhead.
// ═══════════════════════════════════════════════════════════
'use strict';

let _webR    = null;
let _ready   = false;
let _failed  = false;
let _initProm = null;

// ── Initialisation ───────────────────────────────────────────

/**
 * Start loading WebR + flexsurv.
 * @param {(state:'loading'|'installing'|'ready'|'error', msg?:string)=>void} onStatus
 * @returns {Promise<boolean>}  true = ready, false = failed
 */
export function initWebR(onStatus) {
  if (_initProm) return _initProm;

  _initProm = _doInit(onStatus);
  return _initProm;
}

async function _doInit(onStatus) {
  try {
    onStatus('loading');
    const { WebR } = await import('https://webr.r-wasm.org/v0.4.2/webr.mjs');

    _webR = new WebR({ quiet: false });
    await _webR.init();

    onStatus('installing');
    await _webR.installPackages(['flexsurv'], { quiet: true });

    // Pre-load namespace and helper functions once
    await _webR.evalRVoid(`
      suppressPackageStartupMessages({
        library(flexsurv)
        library(survival)
      })

      # Robust gamma fit using MASS::fitdistr equivalent
      .fit_gamma <- function(x) {
        x <- x[x > 0]
        if (length(x) < 3) return(NULL)
        m  <- mean(x)
        v  <- var(x)
        shape0 <- m^2 / v
        rate0  <- m  / v
        tryCatch(
          MASS::fitdistr(x, "gamma",
                         start = list(shape = shape0, rate = rate0),
                         lower = c(1e-4, 1e-6)),
          error = function(e) NULL
        )
      }

      # Single pathway fit: returns named list(p_ks, p_ad)
      .fit_one <- function(null_ks, null_ad, obs_ks, obs_ad, emp_ks, emp_ad) {
        # --- KS: Gamma on |null_ks| ---
        abs_null <- abs(null_ks)
        p_ks <- tryCatch({
          fit <- .fit_gamma(abs_null)
          if (is.null(fit)) {
            emp_ks
          } else {
            pgamma(abs(obs_ks),
                   shape = fit$estimate["shape"],
                   rate  = fit$estimate["rate"],
                   lower.tail = FALSE)
          }
        }, error = function(e) emp_ks)

        # --- AD: Generalised Gamma on scaled null_ad ---
        s_fac  <- mean(null_ad)
        if (s_fac < 1e-12) return(list(p_ks = emp_ks, p_ad = emp_ad))
        s_data <- null_ad / s_fac
        s_obs  <- obs_ad  / s_fac

        p_ad <- tryCatch({
          fit_gg <- flexsurvreg(Surv(s_data) ~ 1, dist = "gengamma")
          mu    <- fit_gg$res["mu",    "est"]
          sigma <- fit_gg$res["sigma", "est"]
          Q     <- fit_gg$res["Q",     "est"]
          pgengamma(s_obs,
                    mu = mu, sigma = sigma, Q = Q,
                    lower.tail = FALSE)
        }, error = function(e) emp_ad)

        list(
          p_ks = max(min(as.numeric(p_ks), 1), 1e-16),
          p_ad = max(min(as.numeric(p_ad), 1), 1e-16)
        )
      }
    `);

    _ready = true;
    onStatus('ready');
    return true;
  } catch (err) {
    _failed = true;
    onStatus('error', String(err.message || err));
    return false;
  }
}

export const isWebRReady  = () => _ready;
export const isWebRFailed = () => _failed;

// ── Per-pathway fitting ──────────────────────────────────────

/**
 * Fit parametric null distributions for ONE pathway.
 *
 * Data is serialized as compact R numeric literals to avoid
 * WebR object-binding API version differences.
 *
 * @param {Float64Array} nullKS
 * @param {Float64Array} nullAD
 * @param {number} obsKS
 * @param {number} obsAD
 * @param {number} empKS  fallback p-value
 * @param {number} empAD  fallback p-value
 * @returns {Promise<{pKS:number, pAD:number, engine:string}>}
 */
export async function fitOne(nullKS, nullAD, obsKS, obsAD, empKS, empAD) {
  if (!_ready) return { pKS: empKS, pAD: empAD, engine: 'empirical' };

  try {
    const rCode = `
      .fit_one(
        null_ks = c(${_vec(nullKS)}),
        null_ad = c(${_vec(nullAD)}),
        obs_ks  = ${_num(obsKS)},
        obs_ad  = ${_num(obsAD)},
        emp_ks  = ${_num(empKS)},
        emp_ad  = ${_num(empAD)}
      )
    `;

    const result = await _webR.evalR(rCode);
    const js     = await result.toJs();

    // toJs() on a named list returns { type:'list', names:[], values:[] }
    const names  = js.names;
    const vals   = js.values;
    const get    = name => {
      const i = names.indexOf(name);
      return i >= 0 ? vals[i].values[0] : null;
    };

    const pKS = get('p_ks');
    const pAD = get('p_ad');

    return {
      pKS: _guardP(pKS, empKS),
      pAD: _guardP(pAD, empAD),
      engine: 'gg+gamma'
    };
  } catch (err) {
    console.warn('[iGSEA] WebR fitOne error:', err);
    return { pKS: empKS, pAD: empAD, engine: 'empirical-fallback' };
  }
}

/**
 * Fit multiple pathways in a single R call.
 * More efficient than looping fitOne() when nPathways > 1.
 *
 * @param {Array<{nullKS, nullAD, obsKS, obsAD, empKS, empAD}>} items
 * @returns {Promise<Array<{pKS, pAD, engine}>>}
 */
export async function fitMany(items) {
  if (!_ready) {
    return items.map(it => ({ pKS: it.empKS, pAD: it.empAD, engine: 'empirical' }));
  }

  try {
    // Build an R list-of-lists call
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
      (function() {
        items <- list(${rItems})
        lapply(items, function(it) {
          .fit_one(it$null_ks, it$null_ad,
                   it$obs_ks,  it$obs_ad,
                   it$emp_ks,  it$emp_ad)
        })
      })()
    `;

    const result = await _webR.evalR(rCode);
    const js     = await result.toJs();

    // js.values is an array of per-pathway list objects
    return js.values.map((pathway, i) => {
      const names = pathway.names;
      const vals  = pathway.values;
      const get   = name => {
        const idx = names.indexOf(name);
        return idx >= 0 ? vals[idx].values[0] : null;
      };
      return {
        pKS:    _guardP(get('p_ks'), items[i].empKS),
        pAD:    _guardP(get('p_ad'), items[i].empAD),
        engine: 'gg+gamma'
      };
    });

  } catch (err) {
    console.warn('[iGSEA] WebR fitMany error:', err);
    return items.map(it => ({ pKS: it.empKS, pAD: it.empAD, engine: 'empirical-fallback' }));
  }
}

// ── Serialisation helpers ────────────────────────────────────

/** Serialize Float64Array to compact R numeric vector string. */
function _vec(arr) {
  // Use fixed notation for values in [0.001, 1e6], else exponential
  const parts = new Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    const v = arr[i];
    parts[i] = (Math.abs(v) < 1e-4 || Math.abs(v) > 1e6)
      ? v.toExponential(6)
      : v.toPrecision(8);
  }
  return parts.join(',');
}

/** Serialize a single number to R literal. */
function _num(v) {
  if (!isFinite(v)) return 'NaN';
  return (Math.abs(v) < 1e-4 || Math.abs(v) > 1e6)
    ? v.toExponential(10)
    : v.toPrecision(12);
}

/** Clamp p-value; fall back to empirical if invalid. */
function _guardP(p, fallback) {
  if (p == null || !isFinite(p) || p < 0 || p > 1) return fallback;
  return Math.max(p, 1e-16);
}
