# 🧬 SNP Analyser

A comprehensive bioinformatics tool for analyzing Single Nucleotide Polymorphisms (SNPs) and their effects on protein sequences. This application provides both command-line and web-based interfaces for DNA sequence analysis, mutation effect prediction, and detailed reporting.

## ✨ Features

- **DNA Sequence Analysis**: Support for both mitochondrial and nuclear DNA sequences
- **Mutation Effect Prediction**: Analyze the impact of SNPs on amino acid sequences
- **Codon Table Support**: Uses appropriate genetic codes (Standard and Vertebrate Mitochondrial)
- **Comprehensive Reporting**: Modular reporting system with multiple analysis types
- **Web Interface**: User-friendly browser-based application
- **CLI Tools**: Command-line interface for programmatic usage
- **Extensible Architecture**: Plugin-based system for adding new analysis modules

## 🚀 Quick Start

### Prerequisites

- Python 3.8+
- BioPython
- FastAPI (for web interface)

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/Sub-department-of-Genetics-ULSL/snp-analyser.git
   cd snp-analyser
   ```

2. **Install dependencies**
   ```bash
   cd backend
   pip install -r requirements.txt
   ```

### Usage Examples

#### Command Line Interface

```python
from analyser_cli import Analyser, Reporter

# Create an analyser instance
analyser = Analyser("ATGCGATCGTAA", dna_type="mitochondrial")

# Define mutations
mutations = [
    {"position": 3, "base": "T"},
    {"position": 8, "base": "A"}
]

# Apply mutations
analyser.apply_mutations(mutations)

# Generate comprehensive report
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

## 📊 Analysis Types

### Mutation Effects Analysis
- **Base Changes**: Track nucleotide substitutions
- **Codon Impact**: Analyze effects on genetic codons
- **Amino Acid Changes**: Predict protein sequence alterations
- **Effect Classification**: 
  - Silent (synonymous) mutations
  - Missense mutations
  - Nonsense (stop codon) mutations
  - Readthrough mutations

### Future Analysis Modules
- **AlphaFold Integration**: Structural impact prediction
- **Conservation Analysis**: Evolutionary conservation scoring
- **Pathogenicity Prediction**: Disease association analysis
- **Population Genetics**: Allele frequency analysis

## 📝 Report Format

The reporting system generates structured, PDF-ready output with consistent formatting:

```
===============================  SNP ANALYSIS REPORT  ===============================

Generated on: 2025-07-23 10:30:15
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

Mutation #1:
  Position: 3
  Change: G → T
  Codon Position: Codon 1, Position 3
  Codon: ATG → ATT
  Amino Acid: M → I
  Effect Type: Missense (amino acid change)

...

------------------------------------------------------------

================================================================================
```

## 🌐 Web Interface Features

- **Interactive Sequence Input**: Paste or upload FASTA sequences
- **Real-time Validation**: Immediate feedback on sequence format
- **Mutation Designer**: Visual mutation specification tool
- **Live Preview**: Real-time analysis results
- **Export Options**: Download reports in multiple formats
- **Responsive Design**: Works on desktop and mobile devices

## 🛠️ Development

### Adding New Analysis Modules

Extend the analysis capabilities by creating new sub-reports:

```python
from analyser_cli import SubReport, ReportFormatter

class MyCustomReport(SubReport):
    @property
    def name(self) -> str:
        return "Custom Analysis"
    
    def generate_content(self, analyser, mutations) -> str:
        # Your analysis logic here
        return ReportFormatter.format_key_value("Result", "Analysis complete")

# Add to reporter
reporter.add_sub_report(MyCustomReport())
```

### Running Tests

```bash
cd backend
python -m pytest tests/
```

## 📄 License

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