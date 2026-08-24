export default function PositionMarker({ genomePosition, genePosition }) {
  return (
    <span className="pos-marker" aria-hidden="true">
      <span className="pos-marker-genome">m.{genomePosition}</span>
      <span className="pos-marker-gene">c.{genePosition}</span>
    </span>
  );
}
