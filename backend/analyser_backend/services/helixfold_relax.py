"""Amber relaxation of a HelixFold-single prediction.

Runs inside the HelixFold-single virtualenv (it needs ``alphafold_paddle`` and
``openmm``), so this module must not import anything from the backend package.
Diagnostics are written to stdout as a single JSON object.
"""

import argparse
import json
import sys

MAX_ITERATIONS = 0
ENERGY_TOLERANCE = 2.39
STIFFNESS = 10.0
MAX_OUTER_ITERATIONS = 3


def main() -> int:
    parser = argparse.ArgumentParser(description="Relax a predicted structure with Amber.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    from alphafold_paddle.common import protein
    from alphafold_paddle.relax import relax

    with open(args.input, encoding="utf-8") as handle:
        prot = protein.from_pdb_string(handle.read())

    relaxer = relax.AmberRelaxation(
        max_iterations=MAX_ITERATIONS,
        tolerance=ENERGY_TOLERANCE,
        stiffness=STIFFNESS,
        exclude_residues=[],
        max_outer_iterations=MAX_OUTER_ITERATIONS,
    )
    relaxed_pdb, debug_data, violations = relaxer.process(prot=prot)

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(relaxed_pdb)

    print(json.dumps({
        "initial_energy": float(debug_data["initial_energy"]),
        "final_energy": float(debug_data["final_energy"]),
        "rmsd": float(debug_data["rmsd"]),
        "attempts": int(debug_data["attempts"]),
        "residue_violations": int(violations.sum()),
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
