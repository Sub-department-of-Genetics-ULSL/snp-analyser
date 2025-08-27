from typing import Dict, List, Optional
from .analyser import Analyser
from .report_base import SubReport, ReportFormatter
from .sub_reporters import MutationEffectsReport


class Reporter:
    """Modular reporter for generating analysis reports from Analyser data."""
    
    def __init__(self):
        """Initialize reporter with default sub-reports."""
        self.sub_reports: Dict[str, SubReport] = {}
        # Add default sub-reports
        self.add_sub_report(MutationEffectsReport())
    
    def add_sub_report(self, sub_report: SubReport) -> None:
        """Add a new sub-report to the reporter.
        
        Args:
            sub_report: SubReport instance to add
        """
        if not isinstance(sub_report, SubReport):
            raise TypeError("sub_report must be an instance of SubReport")
        
        self.sub_reports[sub_report.name] = sub_report
    
    def remove_sub_report(self, name: str) -> None:
        """Remove a sub-report by name.
        
        Args:
            name: Name of the sub-report to remove
        """
        if name in self.sub_reports:
            del self.sub_reports[name]
    
    def list_sub_reports(self) -> List[str]:
        """Return list of available sub-report names."""
        return list(self.sub_reports.keys())
    
    def generate_report(
        self, 
        analyser: Analyser, 
        mutations: List[Dict], 
        include_reports: Optional[List[str]] = None
    ) -> str:
        """Generate complete report with specified sub-reports.
        
        Args:
            analyser: Analyser instance with sequence data
            mutations: List of mutations that were applied
            include_reports: List of sub-report names to include. If None, includes all.
            
        Returns:
            Complete report as a string
        """
        if include_reports is None:
            include_reports = list(self.sub_reports.keys())
        
        # Validate requested reports exist
        missing_reports = set(include_reports) - set(self.sub_reports.keys())
        if missing_reports:
            raise ValueError(f"Unknown sub-reports: {missing_reports}")
        
        report_sections = [
            ReportFormatter.format_section_header("SNP Analysis Report"),
            "",
            ReportFormatter.format_key_value("Generated on", self._get_timestamp()),
            ReportFormatter.format_key_value("Sequence Type", analyser.dna_type.title()),
            ReportFormatter.format_key_value("Original Sequence", str(analyser.input_sequence)),
            "",
            ReportFormatter.format_section_separator(),
            ""
        ]
        
        # Generate each requested sub-report
        for report_name in include_reports:
            if report_name in self.sub_reports:
                sub_report = self.sub_reports[report_name]
                report_sections.extend([
                    sub_report.generate(analyser, mutations),
                    ""
                ])
        
        # Add final separator
        report_sections.append(ReportFormatter.format_section_separator())
        
        return "\n".join(report_sections)
    
    def _get_timestamp(self) -> str:
        """Get current timestamp for report header."""
        from datetime import datetime
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")
