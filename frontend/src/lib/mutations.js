const COMPLEMENT = { A: 'T', C: 'G', G: 'C', T: 'A' };

export function reverseComplement(sequence) {
  return sequence
    .split('')
    .reverse()
    .map((base) => COMPLEMENT[base] || base)
    .join('');
}

/**
 * 1-based genome position (HGVS m.) of a 0-based index in the gene.
 *
 * The gene is shown in transcription order. For a GenBank complement(a..b) location
 * that order runs down the genome: the first base of the gene (c.1) is position b.
 */
export function genomePosition(geneIndex, geneInfo) {
  let position =
    geneInfo.strand === -1 ? geneInfo.endInGenome - geneIndex : geneInfo.startInGenome + geneIndex;
  if (geneInfo.circular && geneInfo.genomeLength) {
    position = ((position - 1 + geneInfo.genomeLength) % geneInfo.genomeLength) + 1;
  }
  return position;
}

/** HGVS c. label of a 0-based gene index: c.-1 before the start, c.*1 after the end. */
export function codingLabel(geneIndex, geneLength) {
  if (geneIndex < 0) {
    return `-${-geneIndex}`;
  }
  if (geneIndex >= geneLength) {
    return `*${geneIndex - geneLength + 1}`;
  }
  return String(geneIndex + 1);
}

export function buildVisualRows(visualRefArray, visualEditArray, geneInfo, rulerInterval = 10) {
  const rows = [];
  let refGeneIndex = 0;
  let newPosition = 1;
  const geneLength = geneInfo.sequenceLength || 0;

  for (let i = 0; i < visualRefArray.length; i++) {
    const refChar = visualRefArray[i];
    const editChar = visualEditArray[i];
    const genome = genomePosition(refGeneIndex, geneInfo);
    const gene = geneInfo.codingNumberingAvailable
      ? codingLabel(refGeneIndex, geneLength)
      : String(refGeneIndex + 1);

    let refGenomeDisplay = '';
    let refGeneDisplay = '';
    if (refChar !== '-') {
      refGenomeDisplay = String(genome);
      refGeneDisplay = gene;
    } else {
      // An inserted base sits in front of the next reference base.
      refGenomeDisplay = `ins before ${genome}`;
      refGeneDisplay = 'ins';
    }

    let newPosDisplay = '';
    if (editChar !== '-') {
      newPosDisplay = String(newPosition);
      newPosition += 1;
    } else {
      newPosDisplay = 'Deleted';
    }

    let refClass = 'nucleotide';
    let editClass = 'nucleotide';
    if (refChar === '-') {
      refClass += ' gap';
      editClass += ' insertion';
    } else if (editChar === '-') {
      editClass += ' deletion';
    } else if (refChar !== editChar) {
      editClass += ' mutated';
    }

    const isRulerRow = refChar !== '-' && rulerInterval > 0 && (refGeneIndex + 1) % rulerInterval === 0;
    if (isRulerRow) {
      refClass += ' ruler-row';
      editClass += ' ruler-row';
    }

    rows.push({
      index: i,
      refChar,
      editChar,
      refClass,
      editClass,
      refGenomeDisplay,
      refGeneDisplay,
      newPosDisplay,
      isRulerRow,
      rulerGenome: isRulerRow ? String(genome) : '',
      rulerGene: isRulerRow ? gene : ''
    });

    if (refChar !== '-') {
      refGeneIndex += 1;
    }
  }

  return rows;
}

/**
 * Group the changed columns of the editor into edit blocks.
 *
 * Changes touching each other form one block, which HGVS describes as a single
 * deletion-insertion. Positions are 0-based indices of the reference gene: the block
 * replaces gene[refStart:refEnd] (empty for an insertion) with `alt`.
 */
export function computeEditBlocks(visualRefArray, visualEditArray) {
  const blocks = [];
  let refIndex = 0;
  let current = null;

  for (let i = 0; i < visualRefArray.length; i++) {
    const refChar = visualRefArray[i];
    const editChar = visualEditArray[i];
    const changed = refChar === '-' || editChar === '-' || refChar !== editChar;

    if (changed) {
      if (!current) {
        current = { firstColumn: i, lastColumn: i, refStart: refIndex, refEnd: refIndex, ref: '', alt: '' };
      }
      current.lastColumn = i;
      if (refChar !== '-') {
        current.ref += refChar;
        current.refEnd = refIndex + 1;
      }
      if (editChar !== '-') {
        current.alt += editChar;
      }
    } else if (current) {
      blocks.push(current);
      current = null;
    }

    if (refChar !== '-') {
      refIndex += 1;
    }
  }

  if (current) {
    blocks.push(current);
  }

  return blocks.filter((block) => block.ref !== block.alt);
}

/**
 * A syntactically valid HGVS description of an edit block as the user made it.
 *
 * It is not normalised: the backend shifts it to the 3' position, recognises
 * duplications, repeats and inversions and merges neighbouring changes.
 */
export function rawDescription(block, geneInfo) {
  const { refStart, refEnd, ref, alt } = block;
  const refLength = refEnd - refStart;

  if (geneInfo.codingNumberingAvailable) {
    const label = (index) => codingLabel(index, geneInfo.sequenceLength);
    if (refLength === 0) {
      return `c.${label(refStart - 1)}_${label(refStart)}ins${alt}`;
    }
    const range = refLength === 1 ? label(refStart) : `${label(refStart)}_${label(refEnd - 1)}`;
    if (!alt) {
      return `c.${range}del`;
    }
    if (refLength === 1 && alt.length === 1) {
      return `c.${label(refStart)}${ref}>${alt}`;
    }
    return `c.${range}delins${alt}`;
  }

  // m. descriptions are written on the plus strand of the genome.
  const minus = geneInfo.strand === -1;
  const plusRef = minus ? reverseComplement(ref) : ref;
  const plusAlt = minus ? reverseComplement(alt) : alt;
  if (refLength === 0) {
    const flanks = [genomePosition(refStart - 1, geneInfo), genomePosition(refStart, geneInfo)];
    const [first, second] = minus ? [flanks[1], flanks[0]] : flanks;
    return `m.${first}_${second}ins${plusAlt}`;
  }

  let first = genomePosition(refStart, geneInfo);
  let last = genomePosition(refEnd - 1, geneInfo);
  if (minus) {
    [first, last] = [last, first];
  }
  const range = refLength === 1 ? `${first}` : `${first}_${last}`;
  if (!alt) {
    return `m.${range}del`;
  }
  if (refLength === 1 && alt.length === 1) {
    return `m.${first}${plusRef}>${plusAlt}`;
  }
  return `m.${range}delins${plusAlt}`;
}

/** Undo edit blocks in the editor arrays. */
export function revertBlocks(visualRefArray, visualEditArray, blocks) {
  const nextRef = [...visualRefArray];
  const nextEdit = [...visualEditArray];
  const ordered = [...blocks].sort((a, b) => b.firstColumn - a.firstColumn);

  for (const block of ordered) {
    for (let column = block.lastColumn; column >= block.firstColumn; column -= 1) {
      if (nextRef[column] === '-') {
        nextRef.splice(column, 1);
        nextEdit.splice(column, 1);
      } else {
        nextEdit[column] = nextRef[column];
      }
    }
  }

  return { visualRefArray: nextRef, visualEditArray: nextEdit };
}

/**
 * Show normalised variants (as returned by the backend) in fresh editor arrays.
 *
 * A variant replaces gene[cdsStart:cdsEnd] with `alt` (coding strand). Replaced bases
 * are shown as substitutions, the rest as deletions or inserted columns.
 */
export function applyVariantsToVisual(baseSequence, variants) {
  const nextRef = baseSequence.split('');
  const nextEdit = baseSequence.split('');
  const ordered = [...variants].sort((a, b) => b.cdsStart - a.cdsStart || b.cdsEnd - a.cdsEnd);

  for (const variant of ordered) {
    const refLength = variant.cdsEnd - variant.cdsStart;
    const shared = Math.min(refLength, variant.alt.length);
    for (let k = 0; k < shared; k += 1) {
      nextEdit[variant.cdsStart + k] = variant.alt[k];
    }
    for (let k = shared; k < refLength; k += 1) {
      nextEdit[variant.cdsStart + k] = '-';
    }

    const inserted = variant.alt.slice(shared).split('');
    if (inserted.length) {
      const column = variant.cdsStart + refLength;
      nextRef.splice(column, 0, ...inserted.map(() => '-'));
      nextEdit.splice(column, 0, ...inserted);
    }
  }

  return { visualRefArray: nextRef, visualEditArray: nextEdit };
}

/**
 * Split pasted text into HGVS descriptions.
 *
 * Whitespace and commas separate descriptions; a semicolon does too, except inside the
 * square brackets of an allele such as m.[8993T>G;9000del].
 */
export function splitDescriptions(text) {
  const tokens = [];
  let current = '';
  let depth = 0;

  for (const character of text) {
    if (character === '[') {
      depth += 1;
    } else if (character === ']') {
      depth = Math.max(depth - 1, 0);
    }

    const separator = /\s/.test(character) || character === ',' || (character === ';' && depth === 0);
    if (separator) {
      if (current) {
        tokens.push(current);
      }
      current = '';
      continue;
    }
    current += character;
  }

  if (current) {
    tokens.push(current);
  }
  return tokens;
}
