// ═══════════════════════════════════════════════════════════
//  core/gsea.js  ·  v2.9 (Original GSEA Empirical FDR)
// ═══════════════════════════════════════════════════════════
'use strict';

import {
  calcSNR, rankOrder, calcGSEAStats,
  cauchyCombine, shuffle,
  calcNES, calcNES_AD, empP_KS, empP_AD
} from './math.js';

import { pvalGG } from './distributions.js';

const MAX_PERMS       = 2000;
const TARGET_CHUNK_MS = 16;

// ── Helper: 二分查找用于快速统计大于等于/小于等于某个值的数量 ──
function countGTE(sortedArr, val) {
  let l = 0, r = sortedArr.length - 1;
  while (l <= r) {
    const m = (l + r) >> 1;
    if (sortedArr[m] >= val) r = m - 1;
    else l = m + 1;
  }
  return sortedArr.length - l;
}

function countLTE(sortedArr, val) {
  let l = 0, r = sortedArr.length - 1;
  while (l <= r) {
    const m = (l + r) >> 1;
    if (sortedArr[m] <= val) l = m + 1;
    else r = m - 1;
  }
  return l;
}

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
      const newChunk = Math.round(chunkSize * TARGET_CHUNK_MS / elapsed);
      chunkSize = Math.max(1, Math.min(500, newChunk));
    }

    const totalElapsed = (performance.now() - t0) / 1000;
    const pct  = Math.round(done / nPerms * 75);
    onProgress(pct, 'Permutations', `${done}/${nPerms} permutations`);

    await _yield();
  }

  if (_ab()) throw new Error('Aborted');

  // ── 3. Empirical p-values & NES ───────────────────────────
  onProgress(76, 'Statistics', 'Empirical p-values & NES…');
  await _yield();

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
  }

  if (_ab()) throw new Error('Aborted');

  // ── 4. AD parametric fitting (engine = 'parametric') ──────
  const useParam    = engine === 'parametric';
  const pKS_arr     = new Array(nP);
  const pAD_arr     = new Array(nP);
  const pAD_par_arr = new Array(nP).fill(null);

  for (let pi = 0; pi < nP; pi++) {
    pKS_arr[pi] = empArr[pi].pKS_emp;
  }

  if (useParam) {
    onProgress(80, 'Fitting', 'Fitting AD null distributions (GG)…');
    await _yield();
    for (let pi = 0; pi < nP; pi++) {
      if (_ab()) throw new Error('Aborted');
      const r   = empArr[pi];
      const res = pvalGG(r.oa, nullAD[pi], r.pAD_emp);
      if (res.fitted) {
        pAD_par_arr[pi] = res.p;
        pAD_arr[pi]     = res.p;
      } else {
        pAD_arr[pi] = r.pAD_emp;
      }
      if (pi % 10 === 0) await _yield();
    }
  } else {
    for (let pi = 0; pi < nP; pi++) pAD_arr[pi] = empArr[pi].pAD_emp;
  }

  // ── 5. Cauchy combination ─────────────────────────────────
  const pCauchy_arr = new Array(nP);
  for (let pi = 0; pi < nP; pi++) {
    pCauchy_arr[pi] = cauchyCombine(pKS_arr[pi], pAD_arr[pi]);
  }

  // ── 6. 提取 GSEA 经验 FDR 需要的所有分布数据 ──────────────
  onProgress(93, 'FDR', 'Building Global NES Distributions…');
  await _yield();

  // 为每个通路计算均值，用于将 Null ES 转化为 Null NES
  const allPermNES_pos = [];
  const allPermNES_neg = [];
  const allPermNES_AD  = [];

  const obsNES_pos = [];
  const obsNES_neg = [];
  const obsNES_AD  = [];

  for (let pi = 0; pi < nP; pi++) {
    let sumPosKS = 0, countPosKS = 0;
    let sumNegKS = 0, countNegKS = 0;
    let sumAD = 0;

    for (let j = 0; j < nPerms; j++) {
      const ks = nullKS[pi][j];
      if (ks >= 0) { sumPosKS += ks; countPosKS++; }
      else         { sumNegKS += ks; countNegKS++; }
      sumAD += nullAD[pi][j];
    }

    const meanPosKS = countPosKS > 0 ? sumPosKS / countPosKS : 0;
    const meanNegKS = countNegKS > 0 ? Math.abs(sumNegKS / countNegKS) : 0;
    const meanAD    = sumAD / nPerms;

    // 构建全排列背景下的 NES
    for (let j = 0; j < nPerms; j++) {
      const ks = nullKS[pi][j];
      if (ks >= 0) { if (meanPosKS > 0) allPermNES_pos.push(ks / meanPosKS); }
      else         { if (meanNegKS > 0) allPermNES_neg.push(ks / meanNegKS); }
      
      if (meanAD > 0) allPermNES_AD.push(nullAD[pi][j] / meanAD);
    }

    // 收集真实观测结果的 NES
    const e = empArr[pi];
    if (e.nes >= 0) obsNES_pos.push(e.nes);
    else obsNES_neg.push(e.nes);
    
    if (e.nes_ad >= 0) obsNES_AD.push(e.nes_ad);
  }

  // 对全分布排序，以便后续进行快速二分查找
  allPermNES_pos.sort((a, b) => a - b);
  allPermNES_neg.sort((a, b) => a - b);
  allPermNES_AD.sort((a, b) => a - b);
  
  obsNES_pos.sort((a, b) => a - b);
  obsNES_neg.sort((a, b) => a - b);
  obsNES_AD.sort((a, b) => a - b);

  // ── 7. 计算原版 GSEA 经验 FDR ───────────────────────────
  onProgress(97, 'FDR', 'Computing Empirical FDRs…');
  await _yield();

  const getFDR = (obsVal, isPositive, obsList, permList) => {
    if (obsList.length === 0 || permList.length === 0) return 1.0;
    
    let probPerm, probObs;
    if (isPositive) {
      probPerm = countGTE(permList, obsVal) / permList.length;
      probObs  = countGTE(obsList, obsVal) / obsList.length;
    } else {
      probPerm = countLTE(permList, obsVal) / permList.length;
      probObs  = countLTE(obsList, obsVal) / obsList.length;
    }
    
    if (probObs === 0) return 1.0;
    return Math.min(1.0, probPerm / probObs); // FDR最高为1
  };

  const results = pathways.map((p, pi) => {
    const e = empArr[pi];
    
    // 计算 KS 对应的经验 FDR (区分正负)
    const fdr_ks = e.nes >= 0 
      ? getFDR(e.nes, true, obsNES_pos, allPermNES_pos)
      : getFDR(e.nes, false, obsNES_neg, allPermNES_neg);

    // 计算 AD 对应的经验 FDR (恒为正向)
    const fdr_ad = getFDR(e.nes_ad, true, obsNES_AD, allPermNES_AD);

    return {
      name:    p.name,
      url:     p.url ?? null,
      size:    p.size,

      es:      e.ok,
      ad:      e.oa,
      nes:     e.nes,
      nes_ad:  e.nes_ad,

      pKS:     pKS_arr[pi],
      pAD:     pAD_arr[pi],
      pCauchy: pCauchy_arr[pi],

      pKS_emp: e.pKS_emp,
      pAD_emp: e.pAD_emp,
      pKS_par: null,
      pAD_par: pAD_par_arr[pi],

      // 分别输出基于 NES 和 NES_AD 计算出的原始机制 FDR
      fdr_ks:  fdr_ks,
      fdr_ad:  fdr_ad,

      curve:   obsStats[pi].curve,
      obsOrd,
      peakIdx: obsStats[pi].peakIdx
    };
  });

  onProgress(100, 'Done', `${results.length} pathway(s) complete`);
  return results;
}

const _yield = () => new Promise(r => setTimeout(r, 0));
