// ═══════════════════════════════════════════════════════════
//  core/gsea.js  ·  v2.9 (Standard GSEA FDR Edition)
// ═══════════════════════════════════════════════════════════
'use strict';

import {
  calcSNR, rankOrder, calcGSEAStats,
  cauchyCombine, shuffle,
  calcNES, calcNES_AD, empP_KS, empP_AD,
  calcGseaFDR // 确保 math.js 中已添加此函数或在本文件末尾补充
} from './math.js';

import { pvalGG } from './distributions.js';

const MAX_PERMS       = 2000;
const TARGET_CHUNK_MS = 10; 

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
    throw new Error(`nCase out of range`);

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

  // ── 2. Permutations ───────────────────────────────────────
  const nullKS = Array.from({ length: nP }, () => new Float64Array(nPerms));
  const nullAD = Array.from({ length: nP }, () => new Float64Array(nPerms));

  const permArr = new Uint16Array(nS);
  for (let i = 0; i < nS; i++) permArr[i] = i;
  const permSNR = new Float32Array(nG);
  const permOrd = new Array(nG);
  const permC   = permArr.subarray(0, nCase);
  const permT   = permArr.subarray(nCase);

  let done = 0;
  let chunkSize = 5;
  const t0 = performance.now();

  while (done < nPerms) {
    if (_ab()) throw new Error('Aborted');
    const tChunk = performance.now();
    const end = Math.min(done + chunkSize, nPerms);

    for (let j = done; j < end; j++) {
      shuffle(permArr);
      calcSNR(exprMat, nG, permC, permT, permSNR);
      rankOrder(permSNR, nG, permOrd);
      for (let pi = 0; pi < nP; pi++) {
        const st = calcGSEAStats(permSNR, permOrd, pathways[pi].mask, nG, false, wt);
        nullKS[pi][j] = st.ks;
        nullAD[pi][j] = st.ad;
      }
    }
    done = end;
    const elapsed = performance.now() - tChunk;
    if (elapsed > 0 && done < nPerms) {
      chunkSize = Math.max(1, Math.min(500, Math.round(chunkSize * TARGET_CHUNK_MS / elapsed)));
    }
    onProgress(Math.round(done/nPerms*70), 'Permutations', `${done}/${nPerms} perms`);
    await _yield();
  }

  // ── 3. Normalization (The GSEA Secret Sauce) ──────────────
  onProgress(71, 'Normalization', 'Calculating Null NES Matrix…');
  
  const nullNES_KS = Array.from({ length: nP }, () => new Float64Array(nPerms));
  const nullNES_AD = Array.from({ length: nP }, () => new Float64Array(nPerms));

  for (let pi = 0; pi < nP; pi++) {
    if (_ab()) throw new Error('Aborted');
    
    // KS 归一化参数：分别计算正负 ES 的均值
    let posSum = 0, posCnt = 0, negSum = 0, negCnt = 0;
    const ksVec = nullKS[pi];
    for (let j = 0; j < nPerms; j++) {
      const v = ksVec[j];
      if (v >= 0) { posSum += v; posCnt++; }
      else { negSum += Math.abs(v); negCnt++; }
    }
    const posMean = posCnt > 0 ? posSum / posCnt : 1e-9;
    const negMean = negCnt > 0 ? negSum / negCnt : 1e-9;

    // AD 归一化参数：计算 AD 的均值 (AD永远为正)
    let adSum = 0;
    const adVec = nullAD[pi];
    for (let j = 0; j < nPerms; j++) adSum += adVec[j];
    const adMean = adSum / nPerms || 1e-9;

    for (let j = 0; j < nPerms; j++) {
      nullNES_KS[pi][j] = ksVec[j] >= 0 ? ksVec[j] / posMean : ksVec[j] / negMean;
      nullNES_AD[pi][j] = adVec[j] / adMean;
    }
    if (pi % 100 === 0) await _yield();
  }

  // ── 4. Empirical Statistics & Observed NES ────────────────
  onProgress(80, 'Statistics', 'p-values & NES…');
  const empArr = new Array(nP);
  for (let pi = 0; pi < nP; pi++) {
    const ok = obsStats[pi].ks;
    const oa = obsStats[pi].ad;
    empArr[pi] = {
      ok, oa,
      pKS_emp: empP_KS(ok, nullKS[pi], nPerms),
      pAD_emp: empP_AD(oa, nullAD[pi], nPerms),
      nes:     calcNES(ok, nullKS[pi]),
      nes_ad:  calcNES_AD(oa, nullAD[pi])
    };
    if (pi % 100 === 99) await _yield();
  }

  // ── 5. Standard GSEA FDR (Pooled Distribution) ────────────
  onProgress(85, 'FDR', 'Pooling distributions…');
  
  const obsNES_KS = empArr.map(e => e.nes);
  const obsNES_AD = empArr.map(e => e.nes_ad);

  // 调用封装在 math.js 中的全局 FDR 计算逻辑
  const fdrKS = internalCalcGseaFDR(obsNES_KS, nullNES_KS, true);
  const fdrAD = internalCalcGseaFDR(obsNES_AD, nullNES_AD, false);

  // ── 6. Parametric Fitting (Optional) ──────────────────────
  const useParam = engine === 'parametric';
  const pAD_final = new Float64Array(nP);
  const pAD_par_arr = new Array(nP).fill(null);

  if (useParam) {
    for (let pi = 0; pi < nP; pi++) {
      if (_ab()) throw new Error('Aborted');
      const res = pvalGG(empArr[pi].oa, nullAD[pi], empArr[pi].pAD_emp);
      if (res.fitted) {
        pAD_par_arr[pi] = res.p;
        pAD_final[pi] = res.p;
      } else {
        pAD_final[pi] = empArr[pi].pAD_emp;
      }
      if (pi % 20 === 0) {
        onProgress(85 + Math.round(pi/nP*10), 'Fitting', `Fitting AD ${pi}/${nP}`);
        await _yield();
      }
    }
  } else {
    for (let pi = 0; pi < nP; pi++) pAD_final[pi] = empArr[pi].pAD_emp;
  }

  // ── 7. Assemble ───────────────────────────────────────────
  onProgress(98, 'Assembling', 'Finalising results…');
  const results = [];
  for (let pi = 0; pi < nP; pi++) {
    const e = empArr[pi];
    const pC = cauchyCombine(e.pKS_emp, pAD_final[pi]);
    
    results.push({
      name:     pathways[pi].name,
      url:      pathways[pi].url ?? null,
      size:     pathways[pi].size,
      es:       e.ok,
      ad:       e.oa,
      nes:      e.nes,
      nes_ad:   e.nes_ad,
      pKS:      e.pKS_emp,
      pAD:      pAD_final[pi],
      pCauchy:  pC,
      fdr_ks:   fdrKS[pi],  // 标准 GSEA FDR (KS)
      fdr_ad:   fdrAD[pi],  // 标准 GSEA FDR (AD)
      pKS_emp:  e.pKS_emp,
      pAD_emp:  e.pAD_emp,
      pAD_par:  pAD_par_arr[pi],
      curve:    obsStats[pi].curve,
      obsOrd,
      peakIdx:  obsStats[pi].peakIdx
    });
  }

  onProgress(100, 'Done', `Complete`);
  return results;
}

/**
 * 内部 GSEA 风格 FDR 计算
 * 如果 math.js 里没写，可以直接放在这里
 */
function internalCalcGseaFDR(obsNES, nullNESMat, twoSided) {
  const nP = obsNES.length;
  const nPerms = nullNESMat[0].length;
  const fdr = new Float64Array(nP);

  if (twoSided) {
    // KS 分正负池
    const nullPos = [], nullNeg = [];
    for (let i = 0; i < nP; i++) {
      const vec = nullNESMat[i];
      for (let j = 0; j < nPerms; j++) {
        if (vec[j] >= 0) nullPos.push(vec[j]); else nullNeg.push(vec[j]);
      }
    }
    const obsPos = obsNES.filter(v => v >= 0);
    const obsNeg = obsNES.filter(v => v < 0);

    for (let i = 0; i < nP; i++) {
      const nes = obsNES[i];
      const isPos = nes >= 0;
      const pool = isPos ? nullPos : nullNeg;
      const obsGroup = isPos ? obsPos : obsNeg;

      const countNull = pool.reduce((acc, v) => isPos ? (v >= nes ? acc + 1 : acc) : (v <= nes ? acc + 1 : acc), 0);
      const countObs  = obsGroup.reduce((acc, v) => isPos ? (v >= nes ? acc + 1 : acc) : (v <= nes ? acc + 1 : acc), 0);
      
      const phiNull = countNull / pool.length;
      const phiObs  = countObs / obsGroup.length;
      fdr[i] = Math.min(1.0, phiNull / (phiObs + 1e-10));
    }
  } else {
    // AD 单向池 (全是正数)
    const nullPool = [];
    for (let i = 0; i < nP; i++) {
      const vec = nullNESMat[i];
      for (let j = 0; j < nPerms; j++) nullPool.push(vec[j]);
    }
    for (let i = 0; i < nP; i++) {
      const nes = obsNES[i];
      const countNull = nullPool.reduce((acc, v) => v >= nes ? acc + 1 : acc, 0);
      const countObs  = obsNES.reduce((acc, v) => v >= nes ? acc + 1 : acc, 0);
      const phiNull = countNull / nullPool.length;
      const phiObs  = countObs / nP;
      fdr[i] = Math.min(1.0, phiNull / (phiObs + 1e-10));
    }
  }
  return fdr;
}

const _yield = () => new Promise(r => setTimeout(r, 0));
