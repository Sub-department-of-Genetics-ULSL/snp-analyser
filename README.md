# 🧬 SNP Analyser

A comprehensive bioinformatics tool for analyzing Single Nucleotide Polymorphisms (SNPs) and their effects on protein sequences. This application provides both command-line and web-based interfaces for DNA sequence analysis, mutation effect prediction, and detailed reporting with support for multiple organisms and their mitochondrial genomes.

## ✨ Features

- **DNA Sequence Analysis**: Support for both mitochondrial and nuclear DNA sequences
- **Multi-Organism Support**: Pre-configured for human, pig, dog, and cow mitochondrial genomes
- **Automatic Genome Management**: Downloads GenBank files from NCBI on demand
- **Mutation Effect Prediction**: Analyze the impact of SNPs on amino acid sequences
- **Multiple Mutation Types**: Support for substitutions, deletions, and insertions
- **Codon Table Support**: Uses appropriate genetic codes (Standard and Vertebrate Mitochondrial)
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

3. **Optional: Install R dependencies for physicochemical analysis**
   ```R
   # In R console
   install.packages("Peptides")
   ```

### Usage Examples

#### Command Line Interface

```python
from backend.analyser_cli import Analyser, Reporter

# Create an analyser instance
analyser = Analyser("ATGCGATCGTAA", dna_type="mitochondrial")

# Define mutations (1-based positions)
mutations = [
    {"type": "sub", "position": 3, "alt": "T"},      # Substitution
    {"type": "del", "position": 8},                  # Deletion
    {"type": "ins", "position": 5, "alt": "ATG"}     # Insertion
]

# Apply mutations
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

# Analyze the gene
analyser = Analyser(sequence, dna_type="mitochondrial")
mutations = [{"type": "sub", "position": 10, "alt": "G"}]
analyser.apply_mutations(mutations)

reporter = Reporter()
report = reporter.generate_report(analyser, mutations)
print(report)
```

#### Web Interface

```bash
# Start the web server
cd backend
uvicorn analyser_backend.main:app --reload

# Open your browser to http://localhost:8000
```

The web interface provides:
- **Organism & Gene Selection**: Choose from pre-configured organisms and their genes
- **Interactive Mutation Designer**: Click on nucleotides to apply mutations
- **Dual-Column Visualization**: Compare reference and mutated sequences side-by-side
- **Coordinate Display**: Shows both gene-relative and genome-absolute positions
- **Real-time Analysis**: Generate reports on-the-fly
- **HGVS-like Notation**: Input mutations using notation like "3G>T", "del.5A", "ins.3_4ATG"


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

## 🌐 Web Interface Features

- **Organism Selection**: Choose from pre-configured organisms (human, pig, dog, cow)
- **Gene Browser**: Select specific genes from mitochondrial genomes
- **Dual-Column Visualization**: Side-by-side comparison of reference and mutated sequences
- **Interactive Mutation Designer**: Click on nucleotides to substitute, delete, or insert
- **Position Tooltips**: Hover to see gene-relative and genome-absolute coordinates
- **HGVS-like Notation Input**: Batch input mutations using notation like "3G>T", "del.5A"
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
- `POST /report/` - Generate analysis report
  ```json
  {
    "organism": "homo_sapiens",
    "gene": "CYTB",
    "mutations": ["14747G>A", "14766C>T", "ins.14750_14751ATG", "del.14755C"]
  }
  ```

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

### Web API Format (HGVS-like)
The web interface accepts HGVS-like notation strings:

- **Substitution**: `"3G>T"` (position 3, G to T)
- **Deletion**: `"del.5A"` or `"5del"` (delete base at position 5)
- **Insertion**: `"ins.10ATG"` or `"ins.10_11ATG"` (insert ATG after position 10)

Positions are specified using **genome-absolute coordinates** and automatically converted to gene-relative positions.

## 🛠️ Development

### Project Structure

```
snp-analyser/
├── backend/
│   ├── analyser_cli/          # Core analysis library
│   │   ├── analyser.py        # Main Analyser class
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

2. The DataManager will automatically download and cache the genome on first use.

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
### DNA Type Selection
- **Mitochondrial**: Uses Vertebrate Mitochondrial codon table (includes UGA as Trp, AGA/AGG as stop)
- **Nuclear**: Uses Standard codon table

Select the appropriate type based on your sequence source:
```python
# For mitochondrial DNA
analyser = Analyser(sequence, dna_type="mitochondrial")

# For nuclear DNA
analyser = Analyser(sequence, dna_type="nuclear")
```

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
analyser = Analyser(gene_data["sequence"], dna_type="mitochondrial")
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