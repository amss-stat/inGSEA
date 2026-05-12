// ═══════════════════════════════════════════════════════════
//  ui/fileio.js  ·  File parsing and drop-zone logic
// ═══════════════════════════════════════════════════════════
'use strict';

/**
 * Parse expression matrix (CSV/TSV/TXT).
 * First row = header (sample names).
 * First col = gene names.
 * Auto log2+1 transforms if max > 20.
 */
export function parseExpr(text) {
  const lines = text.trim().split(/\r?\n/).filter(l => l.trim());
  if (lines.length < 3) throw new Error('Expression file must have ≥3 rows');

  // Detect delimiter
  const tab = (lines[0].match(/\t/g) || []).length;
  const com = (lines[0].match(/,/g)  || []).length;
  const delim = tab >= com ? '\t' : ',';

  const parts0 = lines[0].split(delim).map(cleanCell);
  const sNames = parts0.slice(1);

  const gNames = [], mat = [];

  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split(delim);
    if (parts.length < 2) continue;
    gNames.push(cleanCell(parts[0]));
    const row = new Float64Array(sNames.length);
    for (let j = 0; j < sNames.length; j++) {
      const v = parseFloat(parts[j + 1]);
      row[j] = isNaN(v) ? 0 : v;
    }
    mat.push(row);
  }

  // Auto log2+1
  let mx = 0;
  for (const r of mat) for (const v of r) if (v > mx) mx = v;
  let transformed = false;
  if (mx > 20) {
    for (const r of mat) for (let j = 0; j < r.length; j++) r[j] = Math.log2(r[j] + 1);
    transformed = true;
  }

  return { gNames, sNames, mat, maxRaw: mx, transformed };
}

/**
 * Parse MSigDB GMT format.
 * Tab-delimited: name \t url \t gene1 \t gene2 …
 * URL field is preserved for linking.
 */
export function parseGMT(text) {
  return text.trim().split(/\r?\n/)
    .filter(l => l.trim())
    .map(l => {
      const parts = l.split('\t');
      const name  = cleanCell(parts[0]);
      const url   = cleanCell(parts[1]);   // MSigDB provides URL here
      const genes = parts.slice(2).map(cleanCell).filter(Boolean);
      return { name, url, genes };
    })
    .filter(p => p.genes.length > 0);
}

/**
 * Parse pathway CSV: header row, then Pathway,Genes(;-separated).
 */
export function parsePathwayCSV(text) {
  const lines = text.trim().split(/\r?\n/).filter(l => l.trim());
  const result = [];
  for (let i = 1; i < lines.length; i++) {
    const parts = csvSplit(lines[i]);
    if (parts.length < 2) continue;
    const name  = cleanCell(parts[0]);
    const genes = parts[1].split(';').map(cleanCell).filter(Boolean);
    result.push({ name, url: null, genes });
  }
  return result;
}

/**
 * Build Uint8Array masks for each pathway.
 * Only pathways with ≥5 matched genes are returned.
 */
export function buildMasks(rawPathways, geneNames) {
  const idx = new Map();
  geneNames.forEach((g, i) => idx.set(g, i));
  const nG = geneNames.length;

  return rawPathways.map(p => {
    const mask = new Uint8Array(nG);
    for (const g of p.genes) {
      const i = idx.get(g);
      if (i !== undefined) mask[i] = 1;
    }
    const size = mask.reduce((s, v) => s + v, 0);
    if (size < 5) return null;
    return { name: p.name, url: p.url || null, mask, size };
  }).filter(Boolean);
}

/** Infer MSigDB collection URL from pathway name (fallback heuristic). */
export function inferMSigDBUrl(name, providedUrl) {
  if (providedUrl && providedUrl.startsWith('http')) return providedUrl;
  // MSigDB standard URL pattern
  const base = 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/';
  return base + encodeURIComponent(name) + '.html';
}

// ── Helpers ────────────────────────────────────────────────

function cleanCell(s) {
  if (!s) return '';
  return s.trim().replace(/^["']|["']$/g, '');
}

function csvSplit(line) {
  const out = []; let cur = '', q = false;
  for (const c of line) {
    if (c === '"') { q = !q; }
    else if (c === ',' && !q) { out.push(cur); cur = ''; }
    else cur += c;
  }
  out.push(cur);
  return out;
}

/**
 * Attach drag-and-drop + file-input handlers to a drop zone.
 */
export function setupDropZone(dzId, fiId, onFile) {
  const dz = document.getElementById(dzId);
  const fi = document.getElementById(fiId);

  dz.addEventListener('dragover', e => { e.preventDefault(); dz.classList.add('over'); });
  dz.addEventListener('dragleave', () => dz.classList.remove('over'));
  dz.addEventListener('drop', e => {
    e.preventDefault();
    dz.classList.remove('over');
    const f = e.dataTransfer.files[0];
    if (f) onFile(f);
  });
  fi.addEventListener('change', e => {
    const f = e.target.files[0];
    if (f) onFile(f);
    fi.value = '';   // reset so same file can be re-loaded
  });
}
