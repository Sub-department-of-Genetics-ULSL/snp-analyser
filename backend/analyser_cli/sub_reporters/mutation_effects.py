from typing import Dict, List

from ..analyser import Analyser
from ..hgvs_protein import CdsEdit, ProteinConsequence, predict_consequence, three_letter
from ..report_base import SubReport, ReportFormatter

# CLI mutations use short type names, the web report the HGVS variant types.
TYPE_LABELS = {
    "sub": "substitution",
    "del": "deletion",
    "ins": "insertion",
    "delins": "deletion-insertion",
    "dup": "duplication",
}


class MutationEffectsReport(SubReport):
    """Sub-report for DNA variants and their predicted protein consequences."""

    @property
    def name(self) -> str:
        return "Mutation Effects"

    def generate_content(self, analyser: Analyser, mutations: List[Dict]) -> str:
        """Generate content showing DNA variants and amino acid changes."""
        if not mutations or analyser.mutated_sequence is None:
            return "No mutations applied."

        report_lines = [
            ReportFormatter.format_key_value("DNA Sequence Type", analyser.dna_type.title()),
            ReportFormatter.format_key_value("Genetic Code", analyser.genetic_code_description),
            ReportFormatter.format_key_value("Original Sequence Length", f"{len(analyser.input_sequence)} bp"),
            ReportFormatter.format_key_value("Number of Mutations", str(len(mutations))),
            "",
            ReportFormatter.format_detail_separator(),
            "MUTATION DETAILS:",
            ReportFormatter.format_detail_separator()
        ]

        for number, mutation in enumerate(mutations, 1):
            report_lines.extend(self._describe_mutation(number, analyser, mutation))

        original = str(analyser.amino_acid_translation)
        mutated = str(analyser.mutated_amino_acid_translation or "")
        mutated_length = self._peptide_length(mutated)
        report_lines.extend([
            "",
            "",
            ReportFormatter.format_detail_separator(),
            "AMINO ACID TRANSLATION SUMMARY:",
            ReportFormatter.format_detail_separator(),
            ReportFormatter.format_key_value("Original", original),
            ReportFormatter.format_key_value("Mutated", mutated, indent=1),
            "",
            # Translation ends at the first stop codon, so a premature one shortens the
            # peptide even though the whole translation is shown above.
            ReportFormatter.format_key_value(
                "Peptide Length",
                f"{self._peptide_length(original)} -> {mutated_length} "
                f"amino acid{'' if mutated_length == 1 else 's'} (up to the first stop codon)",
            ),
        ])

        return "\n".join(report_lines)

    @staticmethod
    def _peptide_length(translation: str) -> int:
        return len(translation.split("*", 1)[0])

    def _describe_mutation(self, number: int, analyser: Analyser, mutation: Dict) -> List[str]:
        mut_type = mutation.get("type", "sub")
        kind = TYPE_LABELS.get(mut_type, mut_type)
        lines = ["", f"Mutation #{number} ({kind.upper()}):"]

        for key, label in (("hgvs_m", "HGVS (genome)"), ("hgvs_c", "HGVS (coding DNA)")):
            if mutation.get(key):
                lines.append(ReportFormatter.format_key_value(label, mutation[key], indent=2))

        lines.append(ReportFormatter.format_key_value("Position", self._position(mutation), indent=2))
        lines.append(
            ReportFormatter.format_key_value("Change", self._change(analyser, mutation, kind), indent=2)
        )

        if kind == "substitution":
            lines.extend(self._codon_effect(analyser, mutation["position"]))
        else:
            lines.append(
                ReportFormatter.format_key_value(
                    "Length Change", self._format_length_change(analyser, mutation), indent=2
                )
            )

        protein, consequence = self._protein_change(analyser, mutation)
        if protein:
            lines.append(ReportFormatter.format_key_value("HGVS (protein)", protein, indent=2))
        lines.append(ReportFormatter.format_key_value("Effect Type", consequence, indent=2))
        return lines

    @staticmethod
    def _position(mutation: Dict) -> str:
        position = mutation["position"]
        end_position = mutation.get("end_position", position)
        if mutation.get("type") in {"ins", "insertion"} or (
            mutation.get("hgvs_m") and mutation.get("ref") == ""
        ):
            # An insertion sits between the anchor and the next base.
            return f"{position}_{position + 1}"
        return str(position) if position == end_position else f"{position}_{end_position}"

    def _change(self, analyser: Analyser, mutation: Dict, kind: str) -> str:
        position = mutation["position"]
        end_position = mutation.get("end_position", position)
        reference = mutation.get("ref")
        if reference is None and kind in {"substitution", "deletion", "deletion-insertion"}:
            reference = str(analyser.input_sequence[position - 1:end_position])
        alt = mutation.get("alt") or ""

        if kind == "substitution":
            return f"{reference} -> {alt}"
        if not alt:
            return f"Deletion of {reference}" if reference else "Deletion"
        if not reference:
            word = "Duplication" if kind == "duplication" else "Insertion"
            return f"{word} of {alt}"
        return f"{reference} -> {alt}"

    def _codon_effect(self, analyser: Analyser, position: int) -> List[str]:
        """Codon and amino acid of a substitution, read from the two translations.

        The translations already start with methionine on an alternative start codon and
        complete a partial 3' stop codon with the poly(A) tail, so they are used instead
        of translating the bare codon.
        """
        offset = analyser.coding_offset
        if position <= offset:
            return [ReportFormatter.format_key_value(
                "Codon Position", "5' of the first complete codon", indent=2
            )]

        codon_index = (position - offset - 1) // 3
        codon_start = offset + codon_index * 3
        lines = [ReportFormatter.format_key_value(
            "Codon Position",
            f"Codon {codon_index + 1}, Position {(position - offset - 1) % 3 + 1}",
            indent=2,
        )]

        original_codon = self._codon(str(analyser.input_sequence), codon_start)
        mutated_codon = self._codon(str(analyser.mutated_sequence), codon_start)
        original = str(analyser.amino_acid_translation)
        mutated = str(analyser.mutated_amino_acid_translation or "")
        if codon_index < len(original) and codon_index < len(mutated):
            lines.extend([
                ReportFormatter.format_key_value(
                    "Codon", f"{original_codon} -> {mutated_codon}", indent=2
                ),
                ReportFormatter.format_key_value(
                    "Amino Acid",
                    f"{three_letter(original[codon_index])} -> {three_letter(mutated[codon_index])}",
                    indent=2,
                ),
            ])
        return lines

    @staticmethod
    def _codon(sequence: str, start: int) -> str:
        """A codon; bases completed by the poly(A) tail of the mRNA are in lower case."""
        codon = sequence[start:start + 3]
        return codon + "a" * (3 - len(codon))

    def _protein_change(self, analyser: Analyser, mutation: Dict) -> tuple[str | None, str]:
        if mutation.get("consequence"):
            return mutation.get("hgvs_p"), mutation["consequence"]

        edit = self._cds_edit(analyser, mutation)
        if edit is None:
            return None, "Unknown"
        try:
            consequence: ProteinConsequence = predict_consequence(
                str(analyser.input_sequence),
                [edit],
                analyser.codon_table,
                codon_start=analyser.codon_start,
                transl_except=analyser.transl_except,
            )
        except ValueError:
            return None, "Unknown"
        return consequence.description(), consequence.label

    @staticmethod
    def _cds_edit(analyser: Analyser, mutation: Dict) -> CdsEdit | None:
        mut_type = mutation.get("type", "sub")
        position = mutation["position"]
        end_position = mutation.get("end_position", position)
        alt = (mutation.get("alt") or "").upper()

        if mut_type == "sub":
            return CdsEdit(position - 1, position, alt)
        if mut_type == "del":
            return CdsEdit(position - 1, end_position, "")
        if mut_type == "ins":
            return CdsEdit(position, position, alt)
        if mut_type == "delins":
            return CdsEdit(position - 1, end_position, alt)
        if mut_type == "dup":
            duplicated = alt or str(analyser.input_sequence[position - 1:end_position])
            return CdsEdit(end_position, end_position, duplicated)
        return None

    def _indel_length_change(self, analyser: Analyser, mutation: Dict) -> int | None:
        if "length_change" in mutation:
            return mutation["length_change"]

        edit = self._cds_edit(analyser, mutation)
        if edit is None:
            return None
        return edit.length_change

    def _format_length_change(self, analyser: Analyser, mutation: Dict) -> str:
        length_change = self._indel_length_change(analyser, mutation)
        if length_change is None:
            return "0"
        sign = "+" if length_change > 0 else ""
        frame = "frameshift" if length_change % 3 else "in frame"
        return f"{sign}{length_change} nt ({frame})"
