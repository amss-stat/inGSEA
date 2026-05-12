// ═══════════════════════════════════════════════════════════
//  core/gsea.js  ·  v2.8 (High Responsiveness Edition)
//
//  Key changes from v2.7:
//  • TARGET_CHUNK_MS reduced to 16ms (matches 60fps frame budget).
//  • Empirical p-value loop is now asynchronous with progress updates.
//  • Parametric fitting yields after EVERY pathway for maximum smoothness.
//  • Removed redundant nullAD_arr memory allocation (passes TypedArray directly).
//  • Assembly and FDR stages now include yield points for large pathway sets.
// ═══════════════════════════════════════════════════════════
'use strict';

import {
  calcSNR, rankOrder, calcGSEAStats,
  cauchyCombine, bhFDR, shuffle,
  calcNES, calcNES_AD, empP_KS, empP_AD
} from './math.js';

import { pvalGG } from './distributions.js';

const MAX_PERMS       = 2000;
const TARGET_CHUNK_MS = 10; // 16ms 确保 UI 线程每帧都有机会刷新

export async function runGSEA(opts) {
  const {
    exprMat, geneNames, pathways,
    engine, onProgress, abortSignal
  } = opts;

  const nPerms = Math.min(opts.nPerms | 0, MAX_PERMS);
  const nG     = geneNames.length;
  const nS     = exprMat[0].length;
  const nCase  = opts.nCase | 0;
  const nP     = pathways.length;
  const wt     = typeof opts.weightP === 'number' ? opts.weightP : 1;

  if (nCase < 2 || nCase > nS - 2)
    throw new Error(
      `nCase=${nCase} out of range [2, ${nS - 2}] (nSamples=${nS})`
    );

  const _ab = () => abortSignal?.aborted === true;

  const cIdx = new Uint16Array(nCase);
  const tIdx = new Uint16Array(nS - nCase);
  for (let i = 0; i < nCase; i++)      cIdx[i] = i;
  for (let i = 0; i < nS - nCase; i++) tIdx[i] = i + nCase;

  // ── 1. Observed ───────────────────────────────────────────
  onProgress(0, 'Observed', 'Computing SNR & enrichment…');

  const obsSNR = new Float64Array(nG);
  calcSNR(exprMat, nG, cIdx, tIdx, obsSNR);

  const obsOrd = new Array(nG);
  rankOrder(obsSNR, nG, obsOrd);

  const obsStats = pathways.map(p =>
    calcGSEAStats(obsSNR, obsOrd, p.mask, nG, true, wt)
  );

  if (_ab()) throw new Error('Aborted');

  // ── 2. Permutations ───────────────────────────────────────
  const nullKS = Array.from({ length: nP }, () => new Float64Array(nPerms));
  const nullAD = Array.from({ length: nP }, () => new Float64Array(nPerms));

  const permArr = new Uint16Array(nS);
  for (let i = 0; i < nS; i++) permArr[i] = i;

  const permSNR = new Float32Array(nG);
  const permOrd = new Array(nG);
  const permC   = permArr.subarray(0, nCase);
  const permT   = permArr.subarray(nCase);

  const t0 = performance.now();
  let chunkSize = 5;
  let done = 0;

  while (done < nPerms) {
    if (_ab()) throw new Error('Aborted');

    const end    = Math.min(done + chunkSize, nPerms);
    const tChunk = performance.now();

    for (let j = done; j < end; j++) {
      shuffle(permArr);
      calcSNR(exprMat, nG, permC, permT, permSNR);
      rankOrder(permSNR, nG, permOrd);
      for (let pi = 0; pi < nP; pi++) {
        const st = calcGSEAStats(
          permSNR, permOrd, pathways[pi].mask, nG, false, wt
        );
        nullKS[pi][j] = st.ks;
        nullAD[pi][j] = st.ad;
      }
    }

    done = end;

    const elapsed = performance.now() - tChunk;
    if (elapsed > 0 && done < nPerms) {
      // 动态调整 chunkSize 使得每次占用 CPU 时间接近 TARGET_CHUNK_MS
      const newChunk = Math.round(chunkSize * TARGET_CHUNK_MS / elapsed);
      chunkSize = Math.max(1, Math.min(500, newChunk));
    }

    const totalElapsed = (performance.now() - t0) / 1000;
    const pct  = Math.round(done / nPerms * 80);
    const eta  = totalElapsed > 0.1
      ? ((nPerms - done) / (done / totalElapsed)).toFixed(0) : '—';
    const rate = (done / (totalElapsed + 0.001)).toFixed(0);
    
    onProgress(pct, 'Permutations',
      `${done}/${nPerms} · ETA ${eta}s · ${rate} perm/s`);

    await _yield();
  }

  if (_ab()) throw new Error('Aborted');

  // ── 3. Empirical p & NES ──────────────────────────────────
  onProgress(81, 'Statistics', 'p-values & NES…');
  await _yield();

  // 异步化处理：即使是简单的 map，在 pathway 很多时也会卡顿
  const empArr = new Array(nP);
  for (let pi = 0; pi < nP; pi++) {
    if (_ab()) throw new Error('Aborted');
    
    const ok = obsStats[pi].ks;
    const oa = obsStats[pi].ad;
    empArr[pi] = {
      ok,
      oa,
      pKS_emp: empP_KS(ok, nullKS[pi], nPerms),
      pAD_emp: empP_AD(oa, nullAD[pi], nPerms),
      nes:     calcNES(ok, nullKS[pi]),
      nes_ad:  calcNES_AD(oa, nullAD[pi])
    };

    // 每 50 条路径释放一次主线程
    if (pi % 50 === 49) {
      onProgress(81, 'Statistics', `Empirical p (${pi + 1}/${nP})…`);
      await _yield();
    }
  }

  // ── 4. AD parametric fitting ──────────────────────────────
  const useParam = engine === 'parametric';

  const pKS_arr = empArr.map(e => e.pKS_emp);
  const pAD_arr = new Array(nP);
  const pAD_par_arr = new Array(nP).fill(null);

  if (useParam) {
    let fitSuccess = 0;
    for (let pi = 0; pi < nP; pi++) {
      if (_ab()) throw new Error('Aborted');

      const r = empArr[pi];
      
      // 优化：直接传入 Float64Array (TypedArray)，移除耗时的 Array.from 转换
      const res = pvalGG(r.oa, nullAD[pi], r.pAD_emp);

      if (res.fitted) {
        pAD_par_arr[pi] = res.p;
        pAD_arr[pi]     = res.p;
        fitSuccess++;
      } else {
        pAD_arr[pi]     = r.pAD_emp;
      }

      // 每一个通路都进行 yield，确保在高密度的 MLE 拟合中界面绝对流畅
      const pct = 85 + Math.round((pi / nP) * 10);
      onProgress(pct, 'Fitting', `Fitting AD (${pi + 1}/${nP}) · Success: ${fitSuccess}`);
      await _yield();
    }

    console.log(`GG fit success: ${fitSuccess}/${nP} pathways`);
  } else {
    for (let pi = 0; pi < nP; pi++) {
      pAD_arr[pi] = empArr[pi].pAD_emp;
    }
  }

  // ── 5. Assemble results ───────────────────────────────────
  onProgress(97, 'Assembling', 'Cauchy & FDR…');
  await _yield();

  // 组装结果：也要防止大数组导致的瞬间卡顿
  const results = [];
  for (let pi = 0; pi < nP; pi++) {
    const e  = empArr[pi];
    const pC = cauchyCombine(pKS_arr[pi], pAD_arr[pi]);
    
    results.push({
      name:     pathways[pi].name,
      url:      pathways[pi].url ?? null,
      size:     pathways[pi].size,
      es:       e.ok,
      ad:       e.oa,
      nes:      e.nes,
      nes_ad:   e.nes_ad,
      pKS_emp:  e.pKS_emp,
      pAD_emp:  e.pAD_emp,
      pKS_par:  null,
      pAD_par:  pAD_par_arr[pi],
      pKS:      pKS_arr[pi],
      pAD:      pAD_arr[pi],
      pCauchy:  pC,
      fdr:      null, 
      curve:    obsStats[pi].curve,
      obsOrd,
      peakIdx:  obsStats[pi].peakIdx
    });
    
    if (pi % 100 === 99) await _yield();
  }

  // ── 6. BH-FDR on pCauchy ─────────────────────────────────
  if (results.length >= 2) {
    await _yield(); // FDR 前最后喘息
    const pv = results.map(r => r.pCauchy);
    const qv = bhFDR(pv);
    results.forEach((r, i) => { r.fdr = qv[i]; });
  }

  onProgress(100, 'Done', `${results.length} pathway(s) complete`);
  return results;
}

const _yield = () => new Promise(r => setTimeout(r, 0));
