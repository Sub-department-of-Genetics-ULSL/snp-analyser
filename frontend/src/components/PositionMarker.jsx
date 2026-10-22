export default function PositionMarker({ genomePosition, genePosition, geneLabel = 'c.' }) {
  return (
    <span className="pos-marker" aria-hidden="true">
      <span className="pos-marker-genome">m.{genomePosition}</span>
      <span className="pos-marker-gene">{geneLabel}{genePosition}</span>
    </span>
  );
}
