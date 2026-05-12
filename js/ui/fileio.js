// ═══════════════════════════════════════════════════════════
//  ui/fileio.js  ·  File parsing + drop-zone setup
// ═══════════════════════════════════════════════════════════
'use strict';

// ── Expression matrix ────────────────────────────────────────
/**
 * Parse a CSV/TSV expression matrix.
 * Row 0 = header (sample names, first cell ignored).
 * Col 0 = gene names.
 * Auto log2+1 if max > 20.
 */
export function parseExpr(text) {
  const lines = text.trim().split(/\r?\n/).filter(l => l.trim());
  if (lines.length < 3)
    throw new Error('Expression file needs ≥ 3 rows (header + ≥ 2 genes)');

  const tab = (lines[0].match(/\t/g) ?? []).length;
  const com = (lines[0].match(/,/g)  ?? []).length;
  const delim = tab >= com ? '\t' : ',';

  const h0     = lines[0].split(delim).map(cleanCell);
  const sNames = h0.slice(1);
  const nS     = sNames.length;
  if (nS < 4)
    throw new Error(`Only ${nS} sample columns found. Check delimiter.`);

  const gNames = [];
  const mat    = [];

  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split(delim);
    if (parts.length < 2) continue;
    gNames.push(cleanCell(parts[0]));
    const row = new Float64Array(nS);
    for (let j = 0; j < nS; j++) {
      const v = parseFloat(parts[j + 1]);
      row[j]  = isNaN(v) ? 0 : v;
    }
    mat.push(row);
  }

  if (mat.length === 0)
    throw new Error('No gene rows parsed. Check file format.');

  // Auto log2+1 transform when data looks like raw counts
  let mx = 0;
  for (const r of mat) for (const v of r) if (v > mx) mx = v;
  let transformed = false;
  if (mx > 20) {
    for (const r of mat)
      for (let j = 0; j < r.length; j++)
        r[j] = Math.log2(r[j] + 1);
    transformed = true;
  }

  return { gNames, sNames, mat, maxRaw: mx, transformed };
}

// ── GMT parsing ──────────────────────────────────────────────
/**
 * Parse MSigDB GMT format.
 * Each line: name \t url \t gene1 \t gene2 …
 * The URL in column 1 is preserved for hyperlinks.
 */
export function parseGMT(text) {
  const result = [];
  for (const line of text.trim().split(/\r?\n/)) {
    const t = line.trim();
    if (!t) continue;
    const parts = t.split('\t');
    if (parts.length < 3) continue;
    const name  = cleanCell(parts[0]);
    const url   = cleanCell(parts[1]);     // MSigDB URL — keep as-is
    const genes = parts.slice(2).map(cleanCell).filter(Boolean);
    if (genes.length > 0) result.push({ name, url, genes });
  }
  if (result.length === 0)
    throw new Error('No valid GMT rows found (expect: name\\turl\\tgene1\\t…)');
  return result;
}

// ── Pathway CSV ──────────────────────────────────────────────
/**
 * Parse simple pathway CSV: header, then Pathway,Genes(;-sep).
 */
export function parsePathwayCSV(text) {
  const lines  = text.trim().split(/\r?\n/).filter(l => l.trim());
  const result = [];
  for (let i = 1; i < lines.length; i++) {
    const parts = _csvSplit(lines[i]);
    if (parts.length < 2) continue;
    const name  = cleanCell(parts[0]);
    const genes = parts[1].split(';').map(cleanCell).filter(Boolean);
    if (genes.length > 0) result.push({ name, url: null, genes });
  }
  return result;
}

// ── Mask builder ─────────────────────────────────────────────
/**
 * Build Uint8Array masks (gene in set = 1).
 * Only returns pathways with ≥ 5 matched genes.
 *
 * @param {Array<{name,url,genes}>} rawPathways
 * @param {string[]}                geneNames
 * @returns {Array<{name,url,mask,size}>}
 */
export function buildMasks(rawPathways, geneNames) {
  const idx = new Map();
  geneNames.forEach((g, i) => idx.set(g, i));
  const nG = geneNames.length;

  const out = [];
  for (const p of rawPathways) {
    const mask = new Uint8Array(nG);
    for (const g of p.genes) {
      const i = idx.get(g);
      if (i !== undefined) mask[i] = 1;
    }
    const size = mask.reduce((s, v) => s + v, 0);
    if (size >= 5) out.push({ name: p.name, url: p.url ?? null, mask, size });
  }
  return out;
}

// ── MSigDB URL ───────────────────────────────────────────────
/**
 * Return a usable MSigDB URL for a pathway.
 * If the GMT file provided a real URL, use it.
 * Otherwise construct the standard MSigDB human geneset URL.
 */
export function msigdbUrl(name, providedUrl) {
  if (providedUrl && /^https?:\/\//.test(providedUrl)) return providedUrl;
  return `https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/${encodeURIComponent(name)}.html`;
}

// ── Drop-zone setup ──────────────────────────────────────────
/**
 * Attach drag-and-drop + click-to-browse behaviour to a drop zone.
 * @param {string}           dzId    element id of the .dz div
 * @param {string}           fiId    element id of the hidden <input type=file>
 * @param {(File)=>void}     onFile
 */
export function setupDropZone(dzId, fiId, onFile) {
  const dz = document.getElementById(dzId);
  const fi = document.getElementById(fiId);

  dz.addEventListener('dragover',  e => { e.preventDefault(); dz.classList.add('over'); });
  dz.addEventListener('dragleave', ()  => dz.classList.remove('over'));
  dz.addEventListener('drop', e => {
    e.preventDefault();
    dz.classList.remove('over');
    const f = e.dataTransfer.files[0];
    if (f) onFile(f);
  });
  fi.addEventListener('change', e => {
    const f = e.target.files[0];
    if (f) { onFile(f); fi.value = ''; }
  });
}

// ── Helpers ──────────────────────────────────────────────────
function cleanCell(s) {
  return (s ?? '').trim().replace(/^["']|["']$/g, '');
}

function _csvSplit(line) {
  const out = []; let cur = '', q = false;
  for (const c of line) {
    if (c === '"') q = !q;
    else if (c === ',' && !q) { out.push(cur); cur = ''; }
    else cur += c;
  }
  out.push(cur);
  return out;
}
