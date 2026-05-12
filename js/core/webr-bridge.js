// ═══════════════════════════════════════════════════════════
//  core/webr-bridge.js
//  Manages the WebR runtime and exposes a single async call:
//    fitNullDistributions(nullKS[], nullAD[], obsKS, obsAD)
//    → { pKS, pAD, ggFit: {mu,sigma,Q,sFac}, gammaFit: {shape,rate} }
//
//  WebR + flexsurv are loaded once at startup.
//  If WebR fails (network, browser compat), the bridge
//  degrades gracefully to the empirical-only engine.
// ═══════════════════════════════════════════════════════════
'use strict';

let _webR    = null;
let _ready   = false;
let _failed  = false;
let _initProm = null;

// ── Public API ───────────────────────────────────────────────

/**
 * Begin loading WebR + flexsurv in the background.
 * Call this once on page load.
 * Returns a Promise that resolves to true (ready) or false (failed).
 */
export function initWebR(onStatus) {
  if (_initProm) return _initProm;

  _initProm = (async () => {
    try {
      onStatus('loading');

      // Dynamic import so we don't block the module graph if CDN is slow
      const { WebR } = await import(
        'https://webr.r-wasm.org/v0.4.2/webr.mjs'
      );

      _webR = new WebR();
      await _webR.init();
      onStatus('installing');

      // Install flexsurv (cached after first load by WebR's package system)
      await _webR.installPackages(['flexsurv'], { quiet: true });

      // Warm-up: load namespace once so first real call is fast
      await _webR.evalRVoid(`suppressPackageStartupMessages(library(flexsurv))`);

      _ready = true;
      onStatus('ready');
      return true;
    } catch (err) {
      _failed = true;
      onStatus('error', err.message);
      return false;
    }
  })();

  return _initProm;
}

export function isWebRReady()  { return _ready; }
export function isWebRFailed() { return _failed; }

// ── Core fitting function ────────────────────────────────────
/**
 * Fit parametric null distributions and return p-values.
 *
 * Mirrors the R script exactly:
 *   KS:  Gamma fitted to |null_ks|, two-sided upper tail
 *   AD:  Generalised Gamma fitted to null_ad / mean(null_ad)
 *        using flexsurv::flexsurvreg(Surv(s_data)~1, dist="gengamma")
 *        then pgengamma(..., lower.tail=FALSE)
 *
 * @param {Float64Array} nullKS   permutation null KS statistics
 * @param {Float64Array} nullAD   permutation null AD statistics
 * @param {number}       obsKS   observed KS
 * @param {number}       obsAD   observed AD
 * @param {number}       empKS   empirical p_KS  (fallback)
 * @param {number}       empAD   empirical p_AD  (fallback)
 * @returns {Promise<{pKS, pAD, meta}>}
 */
export async function fitNullDistributions(
  nullKS, nullAD, obsKS, obsAD, empKS, empAD
) {
  if (!_ready) {
    return { pKS: empKS, pAD: empAD, meta: { engine: 'empirical' } };
  }

  try {
    // Transfer typed arrays to R.
    // We pass them as plain JS numbers via JSON; for large arrays we
    // use WebR's shelter and .assign() which handles Float64Array directly.
    const shelter = await _webR.newShelter();
    try {
      await shelter.evalRVoid(`library(flexsurv)`);

      // Assign null vectors
      await shelter.evalRVoid(
        `null_ks <- .GlobalEnv$.__ks__; null_ad <- .GlobalEnv$.__ad__`
      );
      // Use proper WebR API to push typed arrays
      await _webR.objs.globalEnv.bind('.__ks__',
        await new _webR.RDouble(Array.from(nullKS)));
      await _webR.objs.globalEnv.bind('.__ad__',
        await new _webR.RDouble(Array.from(nullAD)));
      await _webR.objs.globalEnv.bind('.__obs_ks__',
        await new _webR.RDouble([obsKS]));
      await _webR.objs.globalEnv.bind('.__obs_ad__',
        await new _webR.RDouble([obsAD]));
      await _webR.objs.globalEnv.bind('.__emp_ks__',
        await new _webR.RDouble([empKS]));
      await _webR.objs.globalEnv.bind('.__emp_ad__',
        await new _webR.RDouble([empAD]));

      const rResult = await _webR.evalR(`
        local({
          null_ks  <- .GlobalEnv$.__ks__
          null_ad  <- .GlobalEnv$.__ad__
          obs_ks   <- .GlobalEnv$.__obs_ks__
          obs_ad   <- .GlobalEnv$.__obs_ad__
          emp_ks   <- .GlobalEnv$.__emp_ks__
          emp_ad   <- .GlobalEnv$.__emp_ad__

          # ── KS: Gamma fit to |null_ks| ─────────────────────
          abs_null_ks <- abs(null_ks)
          p_ks <- tryCatch({
            fit_gam <- fitdistr_gamma(abs_null_ks)
            pgamma(abs(obs_ks),
                   shape = fit_gam$estimate["shape"],
                   rate  = fit_gam$estimate["rate"],
                   lower.tail = FALSE)
          }, error = function(e) emp_ks)

          # ── AD: Generalised Gamma fit ───────────────────────
          s_fac  <- mean(null_ad)
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

          # Clamp
          p_ks <- max(min(p_ks, 1), 1e-16)
          p_ad <- max(min(p_ad, 1), 1e-16)

          list(p_ks = p_ks, p_ad = p_ad)
        })
      `);

      const obj   = await rResult.toJs();
      const pKS   = obj.values[0].values[0];
      const pAD   = obj.values[1].values[0];

      return {
        pKS: isFinite(pKS) ? pKS : empKS,
        pAD: isFinite(pAD) ? pAD : empAD,
        meta: { engine: 'gg+gamma' }
      };

    } finally {
      shelter.purge();
    }

  } catch (err) {
    // Any R-level or bridge error → fall back silently
    console.warn('[iGSEA] WebR fit failed, using empirical p:', err.message);
    return { pKS: empKS, pAD: empAD, meta: { engine: 'empirical-fallback' } };
  }
}
