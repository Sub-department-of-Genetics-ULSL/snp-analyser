# 🧬 SNP Analyser

A comprehensive bioinformatics tool for analyzing Single Nucleotide Polymorphisms (SNPs) and their effects on protein sequences. This application provides both command-line and web-based interfaces for DNA sequence analysis, mutation effect prediction, and detailed reporting with support for multiple organisms and their mitochondrial genomes.

## ✨ Features

- **DNA Sequence Analysis**: Support for both mitochondrial and nuclear DNA sequences
- **Multi-Organism Support**: Pre-configured for human, dog, house mouse, zebrafish, fruit fly, and roundworm mitochondrial genomes
- **Automatic Genome Management**: Downloads GenBank files from NCBI on demand
- **Mutation Effect Prediction**: Analyze the impact of SNPs on amino acid sequences
- **Multiple Mutation Types**: Support for substitutions, deletions, and insertions
- **Codon Table Support**: Uses appropriate genetic codes (Standard, Vertebrate Mitochondrial and Invertebrate Mitochondrial)
- **Comprehensive Reporting**: Modular reporting system with multiple analysis types
- **Physicochemical Properties**: Analyze changes in peptide properties (requires R and Peptides package)
- **Web Interface**: Interactive browser-based application with dual-column gene visualization
- **CLI Tools**: Python API for programmatic usage
- **Extensible Architecture**: Plugin-based system for adding new analysis modules

## 🚀 Quick Start

### Prerequisites

- **Python 3.10+**
- **BioPython** - for sequence analysis
- **FastAPI** - for web interface
- **R (optional)** - for physicochemical properties analysis
  - Peptides package: `install.packages("Peptides")`
  - rpy2: Python-R interface

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/Sub-department-of-Genetics-ULSL/snp-analyser.git
   cd snp-analyser
   ```

2. **Install Python dependencies**
   ```bash
   cd backend
   pip install -r requirements.txt
   ```

3. **Prepare backend configuration**
   ```bash
   cp .env.example .env
   cd ..
   ./scripts/setup_helixfold_single.sh
   ```

   The helper script clones HelixFold-single into `backend/external/HelixFold-single`, downloads the model weights, and prints the exact `backend/.env` values to enable it.
   It is safe to rerun: if the repo or model already exist, the script skips them.
   Then create a clean Python 3.11 venv inside the HelixFold repo, remove version pins from its `requirements.txt`, add `paddlepaddle`, and install the requirements there.

   Example:
   ```bash
   cd backend/external/HelixFold-single
   python3.11 -m venv .venv_backup
   source .venv_backup/bin/activate
   python -m pip install --upgrade pip
   # edit requirements.txt: remove version pins, add paddlepaddle
   python -m pip install -r requirements.txt
   ```

4. **Optional: Install R dependencies for physicochemical analysis**
   ```R
   # In R console
   install.packages("Peptides")
   ```

### Usage Examples

#### Command Line Interface

```python
from backend.analyser_cli import Analyser, Reporter

# Create an analyser instance
analyser = Analyser("ATGCGATCGTAA", transl_table=2)

# Define mutations (1-based positions)
mutations = [
    {"type": "sub", "position": 3, "alt": "T"},      # Substitution
    {"type": "del", "position": 8},                  # Deletion
    {"type": "ins", "position": 5, "alt": "ATG"}     # Insertion
]

# Apply mutations. Malformed input is rejected before anything is applied:
# an unknown type, a missing key, a position outside the sequence or a base
# that is not A/C/G/T raises ValueError (TypeError for wrong types).
analyser.apply_mutations(mutations)

# Generate comprehensive report
reporter = Reporter()
report = reporter.generate_report(analyser, mutations)
print(report)
```

#### Using Genomic Data Manager

```python
from backend.analyser_cli import DataManager, Analyser, Reporter

# Initialize data manager
manager = DataManager()

# List available organisms
print(manager.organisms)  # {'pig': 'sus_scrofa', 'dog': 'canis_lupus_familiaris', ...}

# Get genes for an organism
genes = manager.list_genes_for_organism("homo_sapiens")
print(genes)  # ['CYTB', 'ND1', 'ND2', ...]

# Get sequence for a specific gene
gene_data = manager.get_gene_for_organism("homo_sapiens", "CYTB")
sequence = gene_data["sequence"]
start_pos = gene_data["startInGenome"]

# Analyze the gene with the genetic code annotated on its GenBank CDS feature
analyser = Analyser(
    sequence,
    transl_table=gene_data["translTable"],
    codon_start=gene_data["codonStart"],
    transl_except=gene_data["translExcept"],
)
mutations = [{"type": "sub", "position": 10, "alt": "G"}]
analyser.apply_mutations(mutations)

reporter = Reporter()
report = reporter.generate_report(analyser, mutations)
print(report)
```

#### Web Interface
On MacOS, if you encounter problems with loading R:
```bash
export R_HOME="$(R RHOME)"
export RPY2_CFFI_MODE=ABI
pip uninstall rpy2 rpy2-rinterface rpy2-robjects
pip install --no-binary :all: "rpy2~=3.6"
```
```bash
# Start the web server
cd backend/analyser_backend
uvicorn main:app --reload

# Open your browser to http://localhost:8000
```

The web interface provides:
- **Organism & Gene Selection**: Choose from pre-configured organisms and their genes
- **Interactive Mutation Designer**: Click on nucleotides to apply mutations
- **Dual-Column Visualization**: Compare reference and mutated sequences side-by-side
- **Coordinate Display**: Shows both gene-relative and genome-absolute positions
- **Asynchronous Report Queue**: Reports are generated in background and tracked in a pending queue
- **HTML Reports**: Completed reports open in a new browser tab
- **3D Protein Viewer**: HTML report renders original protein structure from local PDB files (if available)
- **HGVS DNA Notation**: Input mutations using `m.` (genome-absolute) or `c.` (gene-relative) notation like "m.5367C>T", "c.3_4insATG"


## 📊 Analysis Types

### Mutation Effects Analysis
- **Base Changes**: Track nucleotide substitutions, deletions, and insertions
- **Codon Impact**: Analyze effects on genetic codons
- **Amino Acid Changes**: Predict protein sequence alterations
- **Effect Classification**: 
  - Silent (synonymous) mutations
  - Missense mutations
  - Nonsense (stop codon) mutations
  - Readthrough mutations

### Physicochemical Properties Analysis (Optional)
When R and the Peptides package are installed, the reporter includes:
- **Isoelectric Point (pI)**: Change in peptide isoelectric point
- **Molecular Weight**: Change in peptide molecular weight
- **Charge**: Change in peptide net charge
- **Hydrophobicity**: Change in hydrophobicity (Kyte-Doolittle scale)
- **Instability Index**: Change in peptide stability prediction

## 📝 Report Format

The reporting system generates structured, PDF-ready output with consistent formatting:

```
==============================  SNP ANALYSIS REPORT  ==============================

Generated on: 2026-02-19 14:30:15
Sequence Type: Mitochondrial
Original Sequence: ATGCGATCGTAA

================================================================================

------------------------------ Mutation Effects ------------------------------

DNA Sequence Type: Mitochondrial
Original Sequence Length: 12 bp
Number of Mutations: 2

........................................
MUTATION DETAILS:
........................................

Mutation #1 (SUB):
  Position: 3
  Change: G -> T
  Codon Position: Codon 1, Position 3
  Original Codon: ATG
  Mutated Codon: ATT
  Original Amino Acid: M
  Mutated Amino Acid: I
  Effect Type: Missense (amino acid change)

Mutation #2 (DEL):
  Position: 8
  Change: Deletion of C
  Effect Type: Frameshift mutation

........................................

----------------------------- Physicochemical Properties -----------------------------

Isoelectric Point Original Sequence: 5.52
Isoelectric Point Mutated Sequence: 5.88
Change: ↑
........................................
Molecular Weight Original Sequence: 1347.5
Molecular Weight Mutated Sequence: 1289.2
Change: ↓
...

================================================================================
```

### HTML Report

The web interface renders the same analysis as a self-contained HTML page that adds a
**Protein Translation Comparison** section. The original and the mutated translation are
aligned with a global pairwise alignment (identity scoring, gap penalties) and then printed
one under the other:

```
ORIGINAL    1 MPMAN LLLLI VPILI AMAFL MLTER KILGY MQLRK GPNVV GPYGL LQPFA  50
MUTATED     1 MPMAN LLLLI VPI-I AMAFL MLTER KILGY MQLCK GPNVV GPYGL LQPFA  49
                            |                       |
```

- residues are grouped in blocks of five and wrapped into lines of 60, so nothing has to be
  scrolled horizontally,
- every line is numbered on both sides with its own residue positions, which stay correct
  even when an indel shifts the mutated protein,
- a changed residue is highlighted and marked with `|` on the line below,
- `-` marks a gap, i.e. a residue that was inserted or deleted,
- a "Protein Change" card summarises both lengths and the number of changed/unchanged
  residues, and the raw sequences stay available in a collapsed `Raw protein sequences`
  block for copy-pasting.

## 🌐 Web Interface Features

- **Organism Selection**: Choose from pre-configured organisms (human, pig, dog, cow)
- **Gene Browser**: Select specific genes from mitochondrial genomes
- **Dual-Column Visualization**: Side-by-side comparison of reference and mutated sequences
- **Interactive Mutation Designer**: Click on nucleotides to substitute, delete, or insert
- **Position Tooltips**: Hover to see gene-relative and genome-absolute coordinates
- **HGVS DNA Notation Input**: Batch input mutations using `m.` or `c.` notation like "m.5367C>T", "c.3_4insATG"
- **Real-time Analysis**: Generate comprehensive reports on-the-fly
- **Modern Dark UI**: Professional interface optimized for long sequences
- **Responsive Design**: Works on desktop and mobile devices

## 🔧 API Endpoints

The FastAPI backend provides the following endpoints:

### Organisms
- `GET /organisms/` - List all available organisms
- `GET /organisms/{latin_name}/` - List genes for a specific organism
- `GET /organisms/{latin_name}/{gene}/` - Get sequence data for a specific gene

### Report Generation
- `POST /report/` - Generate text analysis report (legacy synchronous endpoint)
- `POST /report/jobs` - Queue asynchronous HTML report generation
- `GET /report/jobs` - List report generation jobs
- `GET /report/jobs/{job_id}` - Get a single job status
  ```json
  {
    "organism": "homo_sapiens",
    "gene": "CYTB",
    "mutations": ["m.14747G>A", "m.14766C>T", "m.14750_14751insATG", "c.5del"]
  }
  ```

The job list is rebuilt from `backend/analyser_backend/generated_reports/` after restart, so completed reports stay visible.
Report job history is also mirrored in `backend/analyser_backend/report_jobs.sqlite3`, which lets the backend restore organism, gene, mutations, and HelixFold flags after a restart.

### Report and PDB storage
- Generated HTML reports are stored in: `backend/analyser_backend/generated_reports/`
- Report metadata SQLite DB: `backend/analyser_backend/report_jobs.sqlite3`
- Local PDB files are served from: `backend/analyser_backend/pdb_files/`
- Expected PDB path per report: `backend/analyser_backend/pdb_files/{organism}/{gene}.pdb`

## 🧪 Mutation Notation

### CLI/Python API Format
Mutations are specified as dictionaries with the following structure:

```python
# Substitution
{"type": "sub", "position": 3, "alt": "T", "ref": "G"}  # ref is optional

# Deletion
{"type": "del", "position": 5, "ref": "A"}  # ref is optional

# Insertion
{"type": "ins", "position": 10, "alt": "ATG"}
```

### Web API Format (HGVS DNA)
The web interface accepts DNA-level HGVS notation with exactly two coordinate prefixes:

- **`m.`** — absolute position in the mitochondrial genome
- **`c.`** — position relative to the start of the gene (1-based)

Every mutation must carry one of these prefixes; unprefixed and legacy forms
(`5del`, `del.5`, `ins.10_11ATG`) are rejected with HTTP 400.

- **Substitution**: `"m.5367C>T"`, `"c.10A>G"`
- **Deletion**: `"m.5444del"`, `"c.10_12del"`
- **Insertion**: `"m.14750_14751insATG"`, `"c.10_11insATG"`
- **Delins**: `"m.14750_14752delinsT"`, `"c.10delinsGG"`
- **Duplication**: `"m.150_152dup"`, `"c.10dup"`

`m.` positions are converted to gene-relative positions using the gene start in the genome
(`gene_position = m_position - gene_start_in_genome + 1`), so `m.` and `c.` notations that point
to the same nucleotide are fully equivalent. Positions resolving outside the selected gene are
rejected. Generated reports always list mutations in genome-absolute `m.` notation and show both
the gene (`c.`) and genome (`m.`) position in separate columns.

## 🛠️ Development

### Project Structure

```
snp-analyser/
├── backend/
│   ├── analyser_cli/          # Core analysis library
│   │   ├── analyser.py        # Main Analyser class
│   │   ├── cds_annotation.py  # GenBank CDS qualifiers and NCBI genetic codes
│   │   ├── protein_alignment.py # Original vs mutated translation alignment
│   │   ├── reporter.py        # Report generation
│   │   ├── data_manager.py    # Genome data management
│   │   ├── fasta_reader.py    # FASTA file parsing
│   │   ├── report_base.py     # Base classes for reports
│   │   ├── variables.py       # Organism/genome mappings
│   │   └── sub_reporters/     # Analysis modules
│   │       ├── mutation_effects.py
│   │       └── physicochemical_properties.py
│   ├── analyser_backend/      # FastAPI web backend
│   │   ├── main.py           # FastAPI app
│   │   ├── routers/          # API endpoints
│   │   ├── templates/        # HTML report template
│   │   └── models/           # Pydantic models
│   ├── tests/                # Test suite
│   └── requirements.txt      # Python dependencies
├── frontend/
│   └── index.html           # Web interface
└── pyproject.toml           # Build configuration
```

### Adding New Analysis Modules

Extend the analysis capabilities by creating new sub-reports:

```python
from backend.analyser_cli.report_base import SubReport, ReportFormatter
from backend.analyser_cli.analyser import Analyser
from typing import Dict, List

class MyCustomReport(SubReport):
    @property
    def name(self) -> str:
        return "Custom Analysis"
    
    def generate_content(self, analyser: Analyser, mutations: List[Dict]) -> str:
        """Generate custom analysis content."""
        if not mutations or analyser.mutated_sequence is None:
            return "No mutations applied."
        
        # Your analysis logic here
        results = []
        results.append(
            ReportFormatter.format_key_value("Analysis Result", "Complete")
        )
        
        return "\n".join(results)

# Add to reporter
from backend.analyser_cli import Reporter
reporter = Reporter()
reporter.add_sub_report(MyCustomReport())
```

### Adding New Organisms

To add support for additional organisms:

1. **Add NCBI ID mapping** in [backend/analyser_cli/variables.py](backend/analyser_cli/variables.py):
   ```python
   LATIN_TO_NCBI_ID_MAPPING = {
       # ... existing entries ...
       "new_organism": "NC_XXXXXX.X",
   }
   
   COMMON_TO_LATIN_MAPPING = {
       # ... existing entries ...
       "common_name": "new_organism",
   }
   ```

2. The DataManager will automatically download and cache the genome on first use. The
   genetic code, reading frame and codon exceptions are read from the GenBank record
   itself (`/transl_table`, `/codon_start`, `/transl_except`), so nothing else has to be
   declared for a new organism.

### Running Tests

```bash
cd backend
python -m pytest tests/ -v
```

### Code Formatting

This project uses [Black](https://github.com/psf/black) for code formatting:

```bash
cd backend
black .
```

Configuration is in [pyproject.toml](pyproject.toml).

## 📦 Dependencies

### Python Packages
- **biopython** (~1.85) - Sequence analysis and GenBank file handling
- **fastapi[standard]** (~0.115.14) - Web framework and server
- **pydantic** - Data validation (included with FastAPI)
- **rpy2** (optional) - Python-R interface for physicochemical properties

### R Packages (Optional)
- **Peptides** - Physicochemical property calculations

### Development
- **black** (~24.0.0) - Code formatter
- **pytest** - Testing framework (in test-requirements.txt)

## ⚙️ Configuration

### NCBI Entrez Email
The DataManager requires an email for NCBI Entrez API access. Set it in [backend/analyser_cli/variables.py](backend/analyser_cli/variables.py):

```python
ENTREZ_EMAIL = "your.email@example.com"
```

### Genomic Data Storage
Downloaded GenBank files are cached in `backend/genomic_data/` to avoid repeated downloads.

### HelixFold-single predictor

The backend loads `backend/.env` automatically on startup.

Use `./scripts/setup_helixfold_single.sh` to clone HelixFold-single into `backend/external/HelixFold-single` and download `helixfold-single.pdparams`.
Create a clean Python 3.11 venv inside the HelixFold repo, strip version pins from its `requirements.txt`, add `paddlepaddle`, and install the requirements there.

Quick check:
```bash
cd backend/external/HelixFold-single
. .venv_backup/bin/activate
python - <<'PY'
import ml_collections, tensorflow, jax, haiku, tree, paddle
print("HelixFold environment OK")
PY
```

Required variables in `backend/.env`:

- `HELIXFOLD_SINGLE_ENABLED=true`
- `HELIXFOLD_SINGLE_REPO_DIR=/absolute/path/to/backend/external/HelixFold-single`
- `HELIXFOLD_SINGLE_MODEL_PATH=/absolute/path/to/backend/models/helixfold-single.pdparams`

Optional variables:

- `HELIXFOLD_SINGLE_PYTHON_BIN=/absolute/path/to/backend/external/HelixFold-single/.venv/bin/python`
- `HELIXFOLD_SINGLE_SCRIPT_RELPATH=helixfold_single_inference.py`
- `HELIXFOLD_SINGLE_TIMEOUT_SECONDS=7200`

If `HELIXFOLD_SINGLE_PYTHON_BIN` is unset, the backend uses the HelixFold repo's `.venv/bin/python` only; it does not fall back to the server interpreter.

When enabled, the report flow runs HelixFold, stores the predicted structure in the report cache, and renders it in the HTML report.

### Genetic Code Selection

The genetic code is never guessed from the organism. It is read per gene from the GenBank
record of the genome, so every one of the NCBI genetic codes (`transl_table` 1-33) is
supported automatically:

| Qualifier | Meaning | Fallback |
| --- | --- | --- |
| `/transl_table` | NCBI genetic code id, e.g. 2 (Vertebrate Mitochondrial) or 5 (Invertebrate Mitochondrial) | 1 (Standard), which is what GenBank means when the qualifier is absent |
| `/codon_start` | 1-based offset of the first complete codon inside the feature | 1 |
| `/transl_except` | Codons whose amino acid differs from the genetic code, used for 3' partial stop codons completed by polyadenylation | none |

`DataManager.get_gene_for_organism` returns them next to the sequence:

```python
manager = DataManager()
gene = manager.get_gene_for_organism("drosophila_melanogaster", "ND1")

analyser = Analyser(
    gene["sequence"],
    transl_table=gene["translTable"],     # 5 for Drosophila
    codon_start=gene["codonStart"],
    transl_except=gene["translExcept"],
)
print(analyser.genetic_code_description)  # Invertebrate Mitochondrial (NCBI transl_table=5)
```

The translation reproduces the GenBank `/translation` qualifier exactly: an alternative
start codon (`ATT`, `GTG`, `TTG`, ...) is reported as methionine and annotated codon
exceptions are applied. `backend/tests/analyser_cli/test_genbank_translation.py` asserts
this for every coding sequence of every cached genome.

Without an explicit id the standard code (NCBI table 1) is used:

```python
analyser = Analyser(sequence)                     # Standard code
analyser = Analyser(sequence, transl_table=2)     # Vertebrate Mitochondrial
```

`dna_type` is only a label shown in the report, it does not influence the translation.

## ⚠️ Important Notes

### Coordinate Systems
- **CLI/Python API**: Uses **1-based, gene-relative** positions
- **Web Interface Input**: Uses **1-based, genome-absolute** positions (automatically converted)
- **Position 1** refers to the first nucleotide in the sequence

### Mutation Processing
- Mutations are applied in **reverse order** by position to maintain correct indices
- Insertions occur **after** the specified position
- Deletions remove the base **at** the specified position
- Invalid mutations (position < 1 or > sequence length) are skipped

### Translation Behavior
- Only complete codons are translated
- If mutated sequence length is not divisible by 3, remaining bases are truncated

### Limitations
- Currently supports only DNA sequences (not RNA)
- Assumes 5' to 3' direction
- No support for ambiguous nucleotides (N, Y, R, etc.)
- Physicochemical analysis requires complete installation of R and rpy2

## 🔬 Example Workflows

### Analyze a specific gene mutation in humans

```python
from backend.analyser_cli import DataManager, Analyser, Reporter

# Setup
manager = DataManager()
reporter = Reporter()

# Get CYTB gene from human mitochondrial genome
gene_data = manager.get_gene_for_organism("homo_sapiens", "CYTB")

# Analyze a known mutation
analyser = Analyser(
    gene_data["sequence"],
    transl_table=gene_data["translTable"],
    codon_start=gene_data["codonStart"],
    transl_except=gene_data["translExcept"],
)
mutations = [
    {"type": "sub", "position": 100, "alt": "T"},
    {"type": "sub", "position": 250, "alt": "G"}
]
analyser.apply_mutations(mutations)

# Generate report
report = reporter.generate_report(analyser, mutations)
print(report)
```

### Compare multiple mutations

```python
from backend.analyser_cli import Analyser, Reporter

sequence = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
reporter = Reporter()

# Test different mutation scenarios
scenarios = [
    {"name": "Single SNP", "mutations": [{"type": "sub", "position": 5, "alt": "T"}]},
    {"name": "Double SNP", "mutations": [{"type": "sub", "position": 5, "alt": "T"}, 
                                         {"type": "sub", "position": 10, "alt": "C"}]},
    {"name": "With Deletion", "mutations": [{"type": "del", "position": 15}]},
]

for scenario in scenarios:
    print(f"\n{'='*60}\n{scenario['name']}\n{'='*60}")
    analyser = Analyser(sequence, dna_type="mitochondrial")
    analyser.apply_mutations(scenario["mutations"])
    report = reporter.generate_report(analyser, scenario["mutations"])
    print(report)
```

## �📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📞 Contact

- **Repository**: [Sub-department-of-Genetics-ULSL/snp-analyser](https://github.com/Sub-department-of-Genetics-ULSL/snp-analyser)
- **Issues**: [GitHub Issues](https://github.com/Sub-department-of-Genetics-ULSL/snp-analyser/issues)

## 🙏 Acknowledgments

- **BioPython**: For sequence analysis capabilities
- **FastAPI**: For modern web API framework
- **Scientific Community**: For genetic code standards and mutation classification systems

---

*Built with ❤️ for the bioinformatics community*