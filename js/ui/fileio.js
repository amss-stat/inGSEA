// ═══════════════════════════════════════════════════════════
//  ui/fileio.js  ·  File parsing + drop zones
// ═══════════════════════════════════════════════════════════
'use strict';

// ── Expression matrix ────────────────────────────────────────
export function parseExpr(text) {
  const lines = text.trim().split(/\r?\n/).filter(l => l.trim());
  if (lines.length < 3)
    throw new Error('Expression file needs ≥ 3 rows (header + ≥ 2 genes)');

  const tab = (lines[0].match(/\t/g) ?? []).length;
  const com = (lines[0].match(/,/g)  ?? []).length;
  const delim = tab >= com ? '\t' : ',';

  const h0     = lines[0].split(delim).map(_clean);
  const sNames = h0.slice(1);
  const nS     = sNames.length;
  if (nS < 4)
    throw new Error(`Only ${nS} sample columns. Check delimiter (tab or comma).`);

  const gNames = [];
  const mat    = [];

  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split(delim);
    if (parts.length < 2) continue;
    gNames.push(_clean(parts[0]));
    const row = new Float64Array(nS);
    for (let j = 0; j < nS; j++) {
      const v = parseFloat(parts[j + 1]);
      row[j]  = isNaN(v) ? 0 : v;
    }
    mat.push(row);
  }

  if (mat.length === 0)
    throw new Error('No gene rows parsed.');

  // Auto log2+1 if max > 20
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

// ── GMT ──────────────────────────────────────────────────────
/**
 * Parse MSigDB GMT format.
 * Each line: name\turl\tgene1\tgene2\t…
 * Works with both .gmt and .txt extensions.
 * No filename restriction — auto-detected by tab count.
 */
export function parseGMT(text) {
  const result = [];
  for (const line of text.trim().split(/\r?\n/)) {
    const t = line.trim();
    if (!t) continue;
    const parts = t.split('\t');
    if (parts.length < 3) continue;   // must have name + url + ≥1 gene
    const name  = _clean(parts[0]);
    const url   = _clean(parts[1]);
    const genes = parts.slice(2).map(_clean).filter(Boolean);
    if (genes.length > 0) result.push({ name, url, genes });
  }
  if (result.length === 0)
    throw new Error(
      'No valid gene set rows found.\n' +
      'Expected GMT format: name⟨tab⟩url⟨tab⟩gene1⟨tab⟩gene2⟨tab⟩…'
    );
  return result;
}

// ── Mask builder ─────────────────────────────────────────────
/**
 * Build masks. Only keeps pathways with ≥ 10 matched genes.
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
    let size = 0;
    for (let i = 0; i < nG; i++) size += mask[i];
    if (size >= 10) out.push({ name: p.name, url: p.url ?? null, mask, size });
  }
  return out;
}

// ── MSigDB URL ───────────────────────────────────────────────
export function msigdbUrl(name, providedUrl) {
  if (providedUrl && /^https?:\/\//.test(providedUrl)) return providedUrl;
  return 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/' +
         encodeURIComponent(name) + '.html';
}

// ── Drop zone ────────────────────────────────────────────────
export function setupDropZone(dzId, fiId, onFile) {
  const dz = document.getElementById(dzId);
  const fi = document.getElementById(fiId);
  dz.addEventListener('dragover',  e => { e.preventDefault(); dz.classList.add('over'); });
  dz.addEventListener('dragleave', ()  => dz.classList.remove('over'));
  dz.addEventListener('drop', e => {
    e.preventDefault(); dz.classList.remove('over');
    if (e.dataTransfer.files[0]) onFile(e.dataTransfer.files[0]);
  });
  fi.addEventListener('change', e => {
    if (e.target.files[0]) { onFile(e.target.files[0]); fi.value = ''; }
  });
}

// ── Helpers ──────────────────────────────────────────────────
function _clean(s) {
  return (s ?? '').trim().replace(/^["']|["']$/g, '');
}
