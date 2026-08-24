export default function Tooltip({ refGeneDisplay, refGenomeDisplay, newPosDisplay }) {
  return (
    <span className="pos-tooltip">
      <div className="tooltip-row"><span className="tooltip-label">Ref Gene Pos (c.):</span> {refGeneDisplay}</div>
      <div className="tooltip-row"><span className="tooltip-label">Ref Genome (m.):</span> {refGenomeDisplay}</div>
      <div className="tooltip-row highlight-row"><span className="tooltip-label">New Pos:</span> {newPosDisplay}</div>
    </span>
  );
}
