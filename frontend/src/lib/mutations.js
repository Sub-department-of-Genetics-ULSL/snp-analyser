export const GENE_PREFIX = 'c';
export const GENOME_PREFIX = 'm';
export const SUPPORTED_PREFIXES = [GENE_PREFIX, GENOME_PREFIX];

/** Convert an HGVS position to a 0-based index inside the gene sequence. */
export function toRelativeIndex(position, prefix, startGenomePos) {
  if (prefix === GENE_PREFIX) {
    return position - 1;
  }
  return position - startGenomePos;
}

/** Convert an HGVS position to an absolute mitochondrial genome position. */
export function toGenomePosition(position, prefix, startGenomePos) {
  return toRelativeIndex(position, prefix, startGenomePos) + startGenomePos;
}

export function buildVisualRows(visualRefArray, visualEditArray, startGenomePos, rulerInterval = 10) {
  const rows = [];
  let refGenomeCounter = startGenomePos;
  let refGeneCounter = 1;
  let newStrandCounter = startGenomePos;

  for (let i = 0; i < visualRefArray.length; i++) {
    const refChar = visualRefArray[i];
    const editChar = visualEditArray[i];

    let refGenomeDisplay = '';
    let refGeneDisplay = '';
    if (refChar !== '-') {
      refGenomeDisplay = String(refGenomeCounter);
      refGeneDisplay = String(refGeneCounter);
    } else {
      refGenomeDisplay = `${refGenomeCounter} (ins)`;
      refGeneDisplay = 'ins';
    }

    let newPosDisplay = '';
    if (editChar !== '-') {
      newPosDisplay = String(newStrandCounter);
      newStrandCounter += 1;
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

    const isRulerRow = refChar !== '-' && rulerInterval > 0 && refGeneCounter % rulerInterval === 0;
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
      rulerGenome: isRulerRow ? String(refGenomeCounter) : '',
      rulerGene: isRulerRow ? String(refGeneCounter) : ''
    });

    if (refChar !== '-') {
      refGenomeCounter += 1;
      refGeneCounter += 1;
    }
  }

  return rows;
}

/** Detected mutations are always reported with genome-absolute `m.` coordinates. */
export function computeMutations(visualRefArray, visualEditArray, startGenomePos) {
  const substitutions = [];
  const deletionPositions = [];
  const insertions = [];
  let currentGenomeIndex = startGenomePos;
  let currentInsertion = null;

  for (let i = 0; i < visualRefArray.length; i++) {
    const refChar = visualRefArray[i];
    const editChar = visualEditArray[i];

    if (refChar === '-') {
      if (!currentInsertion) {
        currentInsertion = {
          left: currentGenomeIndex - 1,
          seq: ''
        };
      }
      currentInsertion.seq += editChar;
    } else {
      if (currentInsertion) {
        insertions.push(currentInsertion);
        currentInsertion = null;
      }

      if (editChar === '-') {
        deletionPositions.push(currentGenomeIndex);
      } else if (refChar !== editChar) {
        substitutions.push(`m.${currentGenomeIndex}${refChar}>${editChar}`);
      }

      currentGenomeIndex += 1;
    }
  }

  if (currentInsertion) {
    insertions.push(currentInsertion);
  }

  const sortedDeletions = [...deletionPositions].sort((a, b) => a - b);
  const deletionBlocks = [];
  for (const position of sortedDeletions) {
    const lastBlock = deletionBlocks[deletionBlocks.length - 1];
    if (lastBlock && position === lastBlock.end + 1) {
      lastBlock.end = position;
    } else {
      deletionBlocks.push({ start: position, end: position });
    }
  }

  const remainingInsertions = [...insertions];

  const mutations = [];

  for (const block of deletionBlocks) {
    const insertionCandidates = remainingInsertions
      .filter((insertion) => insertion.left >= block.start - 1 && insertion.left <= block.end)
      .sort((a, b) => {
        const anchorDistanceA = Math.abs(block.end - a.left);
        const anchorDistanceB = Math.abs(block.end - b.left);
        if (anchorDistanceA !== anchorDistanceB) {
          return anchorDistanceA - anchorDistanceB;
        }
        return b.left - a.left;
      });
    const matchedInsertion = insertionCandidates[0];
    if (matchedInsertion) {
      const insertionIdx = remainingInsertions.indexOf(matchedInsertion);
      if (insertionIdx >= 0) {
        remainingInsertions.splice(insertionIdx, 1);
      }
    }

    const delRange =
      block.start === block.end ? `${block.start}` : `${block.start}_${block.end}`;

    if (matchedInsertion) {
      mutations.push(`m.${delRange}delins${matchedInsertion.seq}`);
    } else {
      mutations.push(`m.${delRange}del`);
    }
  }

  for (const insertion of remainingInsertions) {
    mutations.push(`m.${insertion.left}_${insertion.left + 1}ins${insertion.seq}`);
  }

  mutations.push(...substitutions);

  return mutations.sort((a, b) => {
    const posA = parseInt(a.match(/(\d+)/)?.[1] || '0', 10);
    const posB = parseInt(b.match(/(\d+)/)?.[1] || '0', 10);
    return posA - posB;
  });
}

export function applyTextMutations(baseSequence, textInputMutations, startGenomePos) {
  const nextRef = baseSequence.split('');
  const nextEdit = baseSequence.split('');
  const operations = [];
  const errors = [];

  const substitutionPattern = /^(?<prefix>[cm])\.(?<pos>\d+)(?<ref>[ACGT])>(?<alt>[ACGT])$/i;
  const insertionPattern = /^(?<prefix>[cm])\.(?<left>\d+)_(?<right>\d+)ins(?<seq>[ACGT]+)$/i;
  const deletionPattern = /^(?<prefix>[cm])\.(?<start>\d+)(?:_(?<end>\d+))?del(?:[ACGT]+)?$/i;
  const delinsPattern = /^(?<prefix>[cm])\.(?<start>\d+)(?:_(?<end>\d+))?delins(?<seq>[ACGT]+)$/i;
  const duplicationPattern = /^(?<prefix>[cm])\.(?<start>\d+)(?:_(?<end>\d+))?dup(?:[ACGT]+)?$/i;

  const tokens = textInputMutations.split(/[\s,;]+/).filter(Boolean);
  for (const token of tokens) {
    const subMatch = token.match(substitutionPattern);
    if (subMatch?.groups) {
      operations.push({
        type: 'sub',
        raw: token,
        prefix: subMatch.groups.prefix.toLowerCase(),
        pos: Number(subMatch.groups.pos),
        ref: subMatch.groups.ref.toUpperCase(),
        alt: subMatch.groups.alt.toUpperCase()
      });
      continue;
    }

    const delinsMatch = token.match(delinsPattern);
    if (delinsMatch?.groups) {
      operations.push({
        type: 'delins',
        raw: token,
        prefix: delinsMatch.groups.prefix.toLowerCase(),
        start: Number(delinsMatch.groups.start),
        end: Number(delinsMatch.groups.end || delinsMatch.groups.start),
        seq: delinsMatch.groups.seq.toUpperCase()
      });
      continue;
    }

    const insMatch = token.match(insertionPattern);
    if (insMatch?.groups) {
      operations.push({
        type: 'ins',
        raw: token,
        prefix: insMatch.groups.prefix.toLowerCase(),
        pos: Number(insMatch.groups.left),
        seq: insMatch.groups.seq.toUpperCase()
      });
      continue;
    }

    const delMatch = token.match(deletionPattern);
    if (delMatch?.groups) {
      operations.push({
        type: 'del',
        raw: token,
        prefix: delMatch.groups.prefix.toLowerCase(),
        start: Number(delMatch.groups.start),
        end: Number(delMatch.groups.end || delMatch.groups.start)
      });
      continue;
    }

    const dupMatch = token.match(duplicationPattern);
    if (dupMatch?.groups) {
      operations.push({
        type: 'dup',
        raw: token,
        prefix: dupMatch.groups.prefix.toLowerCase(),
        start: Number(dupMatch.groups.start),
        end: Number(dupMatch.groups.end || dupMatch.groups.start)
      });
      continue;
    }

    errors.push(token);
  }

  const operationPosition = (operation) =>
    toRelativeIndex(operation.start ?? operation.pos, operation.prefix, startGenomePos);

  operations.sort((a, b) => operationPosition(b) - operationPosition(a));

  operations.forEach((op) => {
    if (op.type === 'sub') {
      const relativeIndex = toRelativeIndex(op.pos, op.prefix, startGenomePos);
      if (relativeIndex < 0 || relativeIndex >= nextRef.length) {
        errors.push(op.raw);
        return;
      }
      if (nextRef[relativeIndex] === op.ref) {
        nextEdit[relativeIndex] = op.alt;
      }
    } else if (op.type === 'del' || op.type === 'delins') {
      const startIndex = toRelativeIndex(op.start, op.prefix, startGenomePos);
      const endIndex = toRelativeIndex(op.end, op.prefix, startGenomePos);
      if (startIndex < 0 || endIndex >= nextRef.length) {
        errors.push(op.raw);
        return;
      }
      for (let i = Math.max(startIndex, endIndex); i >= Math.min(startIndex, endIndex); i -= 1) {
        nextEdit[i] = '-';
      }

      if (op.type === 'delins') {
        const insertSeq = op.seq.split('');
        const insertIndex = Math.max(startIndex, endIndex) + 1;
        for (let k = 0; k < insertSeq.length; k += 1) {
          nextRef.splice(insertIndex + k, 0, '-');
          nextEdit.splice(insertIndex + k, 0, insertSeq[k]);
        }
      }
    } else if (op.type === 'ins') {
      const relativeIndex = toRelativeIndex(op.pos, op.prefix, startGenomePos);
      if (relativeIndex < 0 || relativeIndex >= nextRef.length) {
        errors.push(op.raw);
        return;
      }
      const insertSeq = op.seq.split('');
      const insertIndex = relativeIndex + 1;

      for (let k = 0; k < insertSeq.length; k++) {
        nextRef.splice(insertIndex + k, 0, '-');
        nextEdit.splice(insertIndex + k, 0, insertSeq[k]);
      }
    } else if (op.type === 'dup') {
      const startIndex = toRelativeIndex(op.start, op.prefix, startGenomePos);
      const endIndex = toRelativeIndex(op.end, op.prefix, startGenomePos);
      if (startIndex < 0 || endIndex >= nextEdit.length) {
        errors.push(op.raw);
        return;
      }
      const left = Math.min(startIndex, endIndex);
      const right = Math.max(startIndex, endIndex);
      const duplicated = nextEdit.slice(left, right + 1);
      const insertIndex = right + 1;
      for (let k = 0; k < duplicated.length; k += 1) {
        nextRef.splice(insertIndex + k, 0, '-');
        nextEdit.splice(insertIndex + k, 0, duplicated[k]);
      }
    }
  });

  return {
    visualRefArray: nextRef,
    visualEditArray: nextEdit,
    errors
  };
}

export function revertMutation(mutString, visualRefArray, visualEditArray, startGenomePos) {
  const delinsMatch = mutString.match(/^([cm])\.(\d+)(?:_(\d+))?delins([ACGT]+)$/i);
  const insertionMatch = mutString.match(/^([cm])\.(\d+)_(\d+)ins([ACGT]+)$/i);
  const substitutionMatch = mutString.match(/^([cm])\.(\d+)[ACGT]>([ACGT])$/i);
  const deletionMatch = mutString.match(/^([cm])\.(\d+)(?:_(\d+))?del$/i);

  let type = null;
  let prefix = null;
  let start = null;
  let end = null;
  if (delinsMatch) {
    type = 'delins';
    prefix = delinsMatch[1].toLowerCase();
    start = Number(delinsMatch[2]);
    end = Number(delinsMatch[3] || delinsMatch[2]);
  } else if (insertionMatch) {
    type = 'ins';
    prefix = insertionMatch[1].toLowerCase();
    start = Number(insertionMatch[2]);
  } else if (substitutionMatch) {
    type = 'sub';
    prefix = substitutionMatch[1].toLowerCase();
    start = Number(substitutionMatch[2]);
  } else if (deletionMatch) {
    type = 'del';
    prefix = deletionMatch[1].toLowerCase();
    start = Number(deletionMatch[2]);
    end = Number(deletionMatch[3] || deletionMatch[2]);
  }

  if (!type || start === null) {
    return {
      visualRefArray,
      visualEditArray
    };
  }

  start = toGenomePosition(start, prefix, startGenomePos);
  if (end !== null) {
    end = toGenomePosition(end, prefix, startGenomePos);
  }

  const targetPos = start;
  const nextRef = [...visualRefArray];
  const nextEdit = [...visualEditArray];

  if (type === 'ins' || type === 'delins') {
    const left = end === null ? targetPos : Math.min(start, end);
    const right = end === null ? targetPos : Math.max(start, end);
    let currentGenomePos = startGenomePos;
    for (let i = 0; i < nextRef.length; i++) {
      const refChar = nextRef[i];
      const insertionAnchor = currentGenomePos - 1;
      const isMatchingInsertion =
        type === 'ins'
          ? insertionAnchor === targetPos
          : insertionAnchor >= left - 1 && insertionAnchor <= right;

      if (refChar === '-' && isMatchingInsertion) {
        while (i < nextRef.length && nextRef[i] === '-') {
          nextRef.splice(i, 1);
          nextEdit.splice(i, 1);
        }
        break;
      }
      if (refChar !== '-') {
        currentGenomePos += 1;
      }
    }
  }

  let currentGenomePos = startGenomePos;
  const left = end === null ? targetPos : Math.min(start, end);
  const right = end === null ? targetPos : Math.max(start, end);
  for (let i = 0; i < nextRef.length; i++) {
    const refChar = nextRef[i];
    if (refChar !== '-') {
      if (type === 'sub' && currentGenomePos === targetPos) {
        nextEdit[i] = nextRef[i];
        break;
      }
      if ((type === 'del' || type === 'delins') && currentGenomePos >= left && currentGenomePos <= right) {
        nextEdit[i] = nextRef[i];
      }
      currentGenomePos += 1;
    }
  }

  return {
    visualRefArray: nextRef,
    visualEditArray: nextEdit
  };
}
