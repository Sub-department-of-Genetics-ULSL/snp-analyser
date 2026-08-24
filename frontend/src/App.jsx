import { useEffect, useMemo, useRef, useState } from 'react';
import NucleotideCell from './components/NucleotideCell';
import PositionMarker from './components/PositionMarker';
import Tooltip from './components/Tooltip';
import { api } from './api';
import {
  applyTextMutations,
  buildVisualRows,
  computeMutations,
  revertMutation
} from './lib/mutations';

const POLL_INTERVAL_MS = 3000;
const RULER_INTERVAL = 10;

export default function App() {
  const [organisms, setOrganisms] = useState({});
  const [animal, setAnimal] = useState('');
  const [genes, setGenes] = useState([]);
  const [gene, setGene] = useState('');

  const [baseSequence, setBaseSequence] = useState('');
  const [startGenomePos, setStartGenomePos] = useState(0);
  const [geneticCode, setGeneticCode] = useState('');
  const [visualRefArray, setVisualRefArray] = useState([]);
  const [visualEditArray, setVisualEditArray] = useState([]);

  const [textInputMutations, setTextInputMutations] = useState('');
  const [predictMutatedStructure, setPredictMutatedStructure] = useState(false);
  const [isGeneratingReport, setIsGeneratingReport] = useState(false);
  const [isLoadingGenes, setIsLoadingGenes] = useState(false);
  const [isLoadingSequence, setIsLoadingSequence] = useState(false);
  const [activePage, setActivePage] = useState('queue');
  const [reportJobs, setReportJobs] = useState([]);
  const [queueError, setQueueError] = useState('');
  const pendingWindowsRef = useRef(new Map());

  const [menuState, setMenuState] = useState({
    visible: false,
    top: 0,
    left: 0,
    index: null
  });

  const refColRef = useRef(null);
  const editColRef = useRef(null);
  const menuRef = useRef(null);
  const isSyncingScrollRef = useRef(false);
  const scrollTopRef = useRef(0);

  useEffect(() => {
    let mounted = true;

    (async () => {
      try {
        const data = await api.getOrganisms();
        if (mounted) {
          setOrganisms(data);
        }
      } catch (error) {
        console.error(error);
      }
    })();

    return () => {
      mounted = false;
    };
  }, []);

  const refreshJobs = async () => {
    try {
      const data = await api.getReportJobs();
      const jobs = data.jobs || [];
      setReportJobs(jobs);
      setQueueError('');

      for (const job of jobs) {
        const pendingWindow = pendingWindowsRef.current.get(job.id);
        if (!pendingWindow || pendingWindow.closed) {
          pendingWindowsRef.current.delete(job.id);
          continue;
        }

        if (job.status === 'completed' && job.report_url) {
          pendingWindow.location.href = api.toReportUrl(job.report_url);
          pendingWindowsRef.current.delete(job.id);
        } else if (job.status === 'failed') {
          pendingWindow.document.write(
            `<pre style="font-family: monospace; padding: 12px;">Report generation failed:\n${job.error || 'Unknown error'}</pre>`
          );
          pendingWindowsRef.current.delete(job.id);
        }
      }
    } catch (error) {
      console.error(error);
      setQueueError('Failed to refresh report queue.');
    }
  };

  useEffect(() => {
    let active = true;
    const tick = async () => {
      if (!active) {
        return;
      }
      await refreshJobs();
    };

    tick();
    const intervalId = window.setInterval(tick, POLL_INTERVAL_MS);
    return () => {
      active = false;
      window.clearInterval(intervalId);
    };
  }, []);

  useEffect(() => {
    if (activePage !== 'editor') {
      return undefined;
    }

    const refCol = refColRef.current;
    const editCol = editColRef.current;
    if (!refCol || !editCol) {
      return undefined;
    }

    // The editor is unmounted while the queue page is shown, so the columns are
    // fresh DOM nodes every time the user comes back: re-bind and restore the offset.
    refCol.scrollTop = scrollTopRef.current;
    editCol.scrollTop = scrollTopRef.current;

    const syncFromRef = () => {
      if (isSyncingScrollRef.current) {
        return;
      }
      isSyncingScrollRef.current = true;
      editCol.scrollTop = refCol.scrollTop;
      scrollTopRef.current = refCol.scrollTop;
      isSyncingScrollRef.current = false;
    };

    const syncFromEdit = () => {
      if (isSyncingScrollRef.current) {
        return;
      }
      isSyncingScrollRef.current = true;
      refCol.scrollTop = editCol.scrollTop;
      scrollTopRef.current = editCol.scrollTop;
      isSyncingScrollRef.current = false;
    };

    refCol.addEventListener('scroll', syncFromRef);
    editCol.addEventListener('scroll', syncFromEdit);

    return () => {
      refCol.removeEventListener('scroll', syncFromRef);
      editCol.removeEventListener('scroll', syncFromEdit);
    };
  }, [activePage, visualRefArray.length]);

  useEffect(() => {
    const handleOutsideClick = (event) => {
      if (!menuState.visible) {
        return;
      }
      if (menuRef.current && menuRef.current.contains(event.target)) {
        return;
      }
      if (editColRef.current && editColRef.current.contains(event.target)) {
        return;
      }
      setMenuState((prev) => ({ ...prev, visible: false, index: null }));
    };

    window.addEventListener('click', handleOutsideClick);
    return () => window.removeEventListener('click', handleOutsideClick);
  }, [menuState.visible]);

  const visualRows = useMemo(
    () => buildVisualRows(visualRefArray, visualEditArray, startGenomePos, RULER_INTERVAL),
    [visualRefArray, visualEditArray, startGenomePos]
  );

  const mutations = useMemo(
    () => computeMutations(visualRefArray, visualEditArray, startGenomePos),
    [visualRefArray, visualEditArray, startGenomePos]
  );

  const hasLoadedSequence = baseSequence.length > 0;

  const loadGenes = async (latinName) => {
    setAnimal(latinName);
    setGene('');
    setGenes([]);
    setBaseSequence('');
    setGeneticCode('');
    setVisualRefArray([]);
    setVisualEditArray([]);

    if (!latinName) {
      return;
    }

    setIsLoadingGenes(true);
    try {
      const data = await api.getGenes(latinName);
      setGenes(data.genes || []);
    } catch (error) {
      console.error(error);
      setGenes([]);
    } finally {
      setIsLoadingGenes(false);
    }
  };

  const loadSequence = async (latinName, selectedGene) => {
    setGene(selectedGene);

    if (!latinName || !selectedGene) {
      setGeneticCode('');
      return;
    }

    setIsLoadingSequence(true);
    try {
      const data = await api.getGeneSequence(latinName, selectedGene);
      const sequence = data.sequence || '';
      const start = Number(data.startInGenome || 0);

      scrollTopRef.current = 0;
      setBaseSequence(sequence);
      setStartGenomePos(start);
      setGeneticCode(
        data.translTableName ? `${data.translTableName} (transl_table=${data.translTable})` : ''
      );
      setVisualRefArray(sequence.split(''));
      setVisualEditArray(sequence.split(''));
    } catch (error) {
      console.error(error);
      setBaseSequence('');
      setGeneticCode('');
      setVisualRefArray([]);
      setVisualEditArray([]);
    } finally {
      setIsLoadingSequence(false);
    }
  };

  const openMenuForIndex = (event, index) => {
    event.stopPropagation();
    const rect = event.currentTarget.getBoundingClientRect();
    setMenuState({
      visible: true,
      top: window.scrollY + rect.top,
      left: window.scrollX + rect.right + 5,
      index
    });
  };

  const applyMenuAction = (action, value) => {
    if (menuState.index === null) {
      return;
    }

    const nextRef = [...visualRefArray];
    const nextEdit = [...visualEditArray];
    const idx = menuState.index;

    if (action === 'sub') {
      nextEdit[idx] = value;
    } else if (action === 'delete') {
      if (nextRef[idx] === '-') {
        nextRef.splice(idx, 1);
        nextEdit.splice(idx, 1);
      } else {
        nextEdit[idx] = '-';
      }
    } else if (action === 'insert') {
      nextRef.splice(idx + 1, 0, '-');
      nextEdit.splice(idx + 1, 0, 'A');
    }

    setVisualRefArray(nextRef);
    setVisualEditArray(nextEdit);
    setMenuState((prev) => ({ ...prev, visible: false, index: null }));
  };

  const handleRevertMutation = (mutString) => {
    const reverted = revertMutation(mutString, visualRefArray, visualEditArray, startGenomePos);
    setVisualRefArray(reverted.visualRefArray);
    setVisualEditArray(reverted.visualEditArray);
  };

  const handleApplyTextMutations = () => {
    if (!textInputMutations || !baseSequence) {
      return;
    }

    const updated = applyTextMutations(baseSequence, textInputMutations, startGenomePos);
    setVisualRefArray(updated.visualRefArray);
    setVisualEditArray(updated.visualEditArray);

    if (updated.errors.length > 0) {
      alert(
        `Skipped ${updated.errors.length} entry/entries: ${updated.errors.join(', ')}\n\n` +
          'Only c. (gene-relative) and m. (genome-absolute) coordinates are accepted, ' +
          'e.g. m.5367C>T or c.10_12del. Positions must fall inside the selected gene.'
      );
    }
  };

  const generateReport = async () => {
    if (!animal || !gene || mutations.length === 0) {
      alert('No mutations.');
      return;
    }

    const pendingWindow = window.open('', '_blank');
    if (pendingWindow) {
      pendingWindow.document.write('<pre style="font-family: monospace; padding: 12px;">Preparing report job...</pre>');
    }

    setIsGeneratingReport(true);
    try {
      const job = await api.createReportJob(animal, gene, mutations, predictMutatedStructure);
      if (pendingWindow) {
        pendingWindow.document.write(
          `<pre style="font-family: monospace; padding: 12px;">Queued report ${job.id}.\nThis tab will open the report when it is ready.</pre>`
        );
        pendingWindowsRef.current.set(job.id, pendingWindow);
      }
      setActivePage('queue');
      await refreshJobs();
    } catch (error) {
      console.error(error);
      if (pendingWindow) {
        pendingWindow.document.write('<pre style="font-family: monospace; padding: 12px;">Failed to queue report.</pre>');
      }
      alert('Failed to queue report.');
    } finally {
      setIsGeneratingReport(false);
    }
  };

  const openReport = (job) => {
    if (!job.report_url) {
      return;
    }
    window.open(api.toReportUrl(job.report_url), '_blank');
  };

  return (
    <>
      <header>
        <h1>SNP Analyser</h1>
      </header>
      <div className="app-nav">
        <button
          className={activePage === 'queue' ? 'active' : ''}
          onClick={() => setActivePage('queue')}
        >
          Pending Report Generation
        </button>
        <button
          className={activePage === 'editor' ? 'active' : ''}
          onClick={() => setActivePage('editor')}
        >
          Mutation Editor
        </button>
      </div>

      {activePage === 'queue' && (
        <div className="queue-page">
          <div className="queue-header">
            <h2>Pending Report Generation</h2>
            <span className="queue-muted">Refreshes automatically every {POLL_INTERVAL_MS / 1000}s</span>
          </div>
          {queueError && <div className="queue-error">{queueError}</div>}
          <div className="queue-table-wrapper">
            <table className="queue-table">
              <thead>
                <tr>
                  <th>Created</th>
                  <th>Organism</th>
                  <th>Gene</th>
                  <th>Status</th>
                  <th>Mutations</th>
                  <th>Action</th>
                </tr>
              </thead>
              <tbody>
                {reportJobs.length === 0 && (
                  <tr>
                    <td colSpan={6}>No report jobs yet.</td>
                  </tr>
                )}
                {reportJobs.map((job) => (
                  <tr key={job.id}>
                    <td>{new Date(job.created_at).toLocaleString()}</td>
                    <td>{job.organism}</td>
                    <td>{job.gene}</td>
                    <td className={`status-${job.status}`}>{job.status}</td>
                    <td>{(job.mutations || []).length}</td>
                    <td>
                      {job.status === 'completed' && job.report_url ? (
                        <button onClick={() => openReport(job)}>Open HTML Report</button>
                      ) : (
                        <span className="queue-muted">
                          {job.status === 'failed' ? (job.error || 'Failed') : 'Waiting'}
                        </span>
                      )}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {activePage === 'editor' && (
        <div className="main-container">
          <div className="gene-display-area">
            <div className="gene-column" id="ref-col" ref={refColRef}>
              <h2>Reference</h2>
              <div id="reference-gene-content">
                {!hasLoadedSequence && <div className="info-text">Select an animal and gene.</div>}
                {visualRows.map((row) => (
                  <NucleotideCell
                    key={`ref-${row.index}`}
                    klass={row.refClass}
                    value={row.refChar}
                    index={row.index}
                    marker={
                      row.isRulerRow ? (
                        <PositionMarker
                          genomePosition={row.rulerGenome}
                          genePosition={row.rulerGene}
                        />
                      ) : null
                    }
                    tooltip={
                      <Tooltip
                        refGeneDisplay={row.refGeneDisplay}
                        refGenomeDisplay={row.refGenomeDisplay}
                        newPosDisplay={row.newPosDisplay}
                      />
                    }
                  />
                ))}
              </div>
            </div>

            <div className="gene-column" id="edit-col" ref={editColRef}>
              <h2>Analysed</h2>
              <div id="editable-gene-content">
                {isLoadingSequence && <div className="info-text">Loading...</div>}
                {!isLoadingSequence &&
                  visualRows.map((row) => (
                    <NucleotideCell
                      key={`edit-${row.index}`}
                      klass={row.editClass}
                      value={row.editChar}
                      index={row.index}
                      onClick={(event) => openMenuForIndex(event, row.index)}
                      marker={
                        row.isRulerRow ? (
                          <PositionMarker
                            genomePosition={row.rulerGenome}
                            genePosition={row.rulerGene}
                          />
                        ) : null
                      }
                      tooltip={
                        <Tooltip
                          refGeneDisplay={row.refGeneDisplay}
                          refGenomeDisplay={row.refGenomeDisplay}
                          newPosDisplay={row.newPosDisplay}
                        />
                      }
                    />
                  ))}
              </div>
            </div>
          </div>

          <div className="controls-panel">
            <div className="control-group">
              <label htmlFor="animal-select">1. Select Animal</label>
              <select id="animal-select" value={animal} onChange={(event) => loadGenes(event.target.value)}>
                <option value="">{Object.keys(organisms).length ? 'Select an animal' : 'Loading...'}</option>
                {Object.entries(organisms).map(([displayName, latinName]) => (
                  <option key={latinName} value={latinName}>
                    {`${displayName.toUpperCase()} ${String(latinName).replaceAll('_', ' ')}`}
                  </option>
                ))}
              </select>
            </div>

            <div className="control-group">
              <label htmlFor="gene-select">2. Select Gene</label>
              <select
                id="gene-select"
                disabled={!animal || isLoadingGenes}
                value={gene}
                onChange={(event) => loadSequence(animal, event.target.value)}
              >
                <option value="">
                  {isLoadingGenes ? 'Loading...' : animal ? 'Select a gene' : 'Select an animal first'}
                </option>
                {genes.map((item) => (
                  <option key={item} value={item}>
                    {item}
                  </option>
                ))}
              </select>
              {geneticCode && (
                <p className="info-text">
                  Genetic code: <strong>{geneticCode}</strong> (from the GenBank{' '}
                  <code>/transl_table</code> qualifier of this gene)
                </p>
              )}
            </div>

            <hr style={{ borderColor: 'var(--border-color)', borderStyle: 'dashed' }} />

            <div className="control-group">
              <label>3. Detected Mutations</label>
              <ul id="mutations-list">
                {mutations.length === 0 && <li>No mutations detected.</li>}
                {mutations.map((mutation) => (
                  <li key={mutation}>
                    <span>{mutation}</span>
                    <button
                      className="remove-mut-btn"
                      title="Revert this mutation"
                      onClick={() => handleRevertMutation(mutation)}
                    >
                      &times;
                    </button>
                  </li>
                ))}
              </ul>
            </div>

            <div className="control-group">
              <label htmlFor="text-input-mutations">Batch Input (Paste text, will remove already present mutations)</label>
              <textarea
                id="text-input-mutations"
                placeholder="e.g. m.5367C>T, c.10_12del, m.21_22insAG, c.100_101delinsTT"
                value={textInputMutations}
                onChange={(event) => setTextInputMutations(event.target.value)}
              />
              <p className="info-text">
                Only <code>m.</code> (absolute position in the mitochondrial genome) and{' '}
                <code>c.</code> (position relative to the gene start) coordinates are accepted.
              </p>
              <button style={{ marginTop: '8px' }} onClick={handleApplyTextMutations}>
                Apply Text Mutations
              </button>
            </div>

            <hr style={{ borderColor: 'var(--border-color)', borderStyle: 'dashed' }} />

            <div className="control-group">
              <label className="checkbox-row" htmlFor="predict-mutated-structure">
                <div>
                  <strong>Predict mutated protein with HelixFold</strong>
                  <p className="info-text">
                    This can take a long time, but cached PDB files will be reused for the same
                    organism, gene, and mutations.
                  </p>
                </div>
                <input
                  checked={predictMutatedStructure}
                  id="predict-mutated-structure"
                  onChange={(event) => setPredictMutatedStructure(event.target.checked)}
                  type="checkbox"
                />
              </label>
            </div>

            <div className="control-group">
              <button disabled={!hasLoadedSequence || isGeneratingReport} onClick={generateReport}>
                {isGeneratingReport ? 'Queuing...' : 'Generate HTML Report'}
              </button>
              <button
                style={{ marginTop: '8px' }}
                className="secondary-btn"
                onClick={() => setActivePage('queue')}
              >
                Open Pending Queue
              </button>
            </div>
          </div>
        </div>
      )}

      {activePage === 'editor' && (
        <div
          id="nucleotide-menu"
          className={menuState.visible ? 'visible' : ''}
          style={{ top: `${menuState.top}px`, left: `${menuState.left}px` }}
          ref={menuRef}
        >
          <div className="menu-row">
            <button data-action="sub" data-value="A" onClick={() => applyMenuAction('sub', 'A')}>
              A
            </button>
            <button data-action="sub" data-value="C" onClick={() => applyMenuAction('sub', 'C')}>
              C
            </button>
          </div>
          <div className="menu-row">
            <button data-action="sub" data-value="G" onClick={() => applyMenuAction('sub', 'G')}>
              G
            </button>
            <button data-action="sub" data-value="T" onClick={() => applyMenuAction('sub', 'T')}>
              T
            </button>
          </div>
          <button id="btn-insert" data-action="insert" onClick={() => applyMenuAction('insert')}>
            Insert Gap
          </button>
          <button id="btn-delete" data-action="delete" onClick={() => applyMenuAction('delete')}>
            Delete
          </button>
        </div>
      )}
    </>
  );
}
