export default function NucleotideCell({ klass, value, tooltip, marker, onClick, index }) {
  return (
    <div className={klass} data-index={index} onClick={onClick}>
      {value}
      {marker}
      {tooltip}
    </div>
  );
}
