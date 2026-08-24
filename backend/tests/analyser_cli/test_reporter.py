"""
Test module for reporter.py classes.

This module provides comprehensive tests for:
- ReportFormatter: Utility class for consistent report formatting
- SubReport: Abstract base class for sub-reports
- Reporter: Main reporter class for generating analysis reports
"""

import pytest
from unittest.mock import Mock, patch
from abc import ABC
from datetime import datetime
from analyser_cli.report_base import ReportFormatter, SubReport
from analyser_cli.reporter import Reporter
from analyser_cli.analyser import Analyser


class TestReportFormatter:
    """Test suite for ReportFormatter utility class."""
    
    def test_format_section_header(self):
        """Test section header formatting with proper padding."""
        header = ReportFormatter.format_section_header("test section")
        
        # Should contain the title in uppercase
        assert "TEST SECTION" in header
        # Should use section characters
        assert ReportFormatter.SECTION_CHAR in header
        # Should be the correct width
        assert len(header) == ReportFormatter.SECTION_WIDTH
    
    def test_format_subsection_header(self):
        """Test subsection header formatting."""
        header = ReportFormatter.format_subsection_header("Test Subsection")
        
        assert "Test Subsection" in header
        assert ReportFormatter.SUBSECTION_CHAR in header
        assert len(header) == ReportFormatter.SUBSECTION_WIDTH
    
    def test_format_section_separator(self):
        """Test section separator line creation."""
        separator = ReportFormatter.format_section_separator()
        
        assert separator == ReportFormatter.SECTION_CHAR * ReportFormatter.SECTION_WIDTH
        assert len(separator) == ReportFormatter.SECTION_WIDTH
    
    def test_format_subsection_separator(self):
        """Test subsection separator line creation."""
        separator = ReportFormatter.format_subsection_separator()
        
        assert separator == ReportFormatter.SUBSECTION_CHAR * ReportFormatter.SUBSECTION_WIDTH
        assert len(separator) == ReportFormatter.SUBSECTION_WIDTH
    
    def test_format_detail_separator(self):
        """Test detail separator line creation."""
        separator = ReportFormatter.format_detail_separator()
        
        assert separator == ReportFormatter.DETAIL_CHAR * ReportFormatter.DETAIL_WIDTH
        assert len(separator) == ReportFormatter.DETAIL_WIDTH
    
    def test_format_key_value_no_indent(self):
        """Test key-value formatting without indentation."""
        result = ReportFormatter.format_key_value("Name", "Value")
        assert result == "Name: Value"
    
    def test_format_key_value_with_indent(self):
        """Test key-value formatting with indentation."""
        result = ReportFormatter.format_key_value("Name", "Value", indent=4)
        assert result == "    Name: Value"
    
    def test_format_list_item_default(self):
        """Test list item formatting with default bullet and indent."""
        result = ReportFormatter.format_list_item("Test item")
        assert result == "  • Test item"
    
    def test_format_list_item_custom(self):
        """Test list item formatting with custom bullet and indent."""
        result = ReportFormatter.format_list_item("Test item", bullet="-", indent=4)
        assert result == "    - Test item"
    
    def test_wrap_content_block_with_title(self):
        """Test content block wrapping with title."""
        content = "Test content line"
        result = ReportFormatter.wrap_content_block(content, "Test Title")
        
        lines = result.split('\n')
        # Should have title header
        assert "Test Title" in lines[0]
        # Should have content
        assert content in lines
        # Should have separators
        assert ReportFormatter.SUBSECTION_CHAR * ReportFormatter.SUBSECTION_WIDTH in lines
    
    def test_wrap_content_block_without_title(self):
        """Test content block wrapping without title."""
        content = "Test content line"
        result = ReportFormatter.wrap_content_block(content)
        
        lines = result.split('\n')
        # Should start with separator
        assert lines[0] == ReportFormatter.SUBSECTION_CHAR * ReportFormatter.SUBSECTION_WIDTH
        # Should contain content
        assert content in lines


class MockSubReport(SubReport):
    """Mock implementation of SubReport for testing."""
    
    def __init__(self, name="Test Report", content="Test content"):
        self._name = name
        self._content = content
    
    @property
    def name(self) -> str:
        return self._name
    
    def generate_content(self, analyser, mutations):
        return self._content


class TestSubReport:
    """Test suite for SubReport abstract base class."""
    
    @pytest.fixture
    def mock_sub_report(self):
        """Create a mock SubReport instance."""
        return MockSubReport()
    
    @pytest.fixture
    def sample_analyser(self):
        """Create a sample Analyser instance for testing."""
        return Analyser("ATGCGATCG", dna_type="mitochondrial")
    
    @pytest.fixture
    def sample_mutations(self):
        """Create sample mutations for testing."""
        return [{"position": 3, "alt": "T"}]
    
    def test_mock_sub_report_properties(self, mock_sub_report):
        """Test mock SubReport properties."""
        assert mock_sub_report.name == "Test Report"
        assert hasattr(mock_sub_report, 'generate_content')
        assert hasattr(mock_sub_report, 'generate')
    
    def test_generate_content(self, mock_sub_report, sample_analyser, sample_mutations):
        """Test generate_content method implementation."""
        content = mock_sub_report.generate_content(sample_analyser, sample_mutations)
        assert content == "Test content"
    
    def test_generate_with_formatting(self, mock_sub_report, sample_analyser, sample_mutations):
        """Test generate method applies proper formatting."""
        result = mock_sub_report.generate(sample_analyser, sample_mutations)
        
        # Should contain the report name
        assert "Test Report" in result
        # Should contain the content
        assert "Test content" in result
        # Should have formatting separators
        assert ReportFormatter.SUBSECTION_CHAR in result


class TestReporter:
    """Test suite for Reporter class."""
    
    @pytest.fixture
    def reporter(self):
        """Create a Reporter instance."""
        return Reporter()
    
    @pytest.fixture
    def sample_analyser(self):
        """Create a sample Analyser instance for testing."""
        analyser = Analyser("ATGCGATCG", dna_type="mitochondrial")
        return analyser
    
    @pytest.fixture
    def sample_mutations(self):
        """Create sample mutations for testing."""
        return [{"position": 3, "alt": "T"}]
    
    @pytest.fixture
    def mock_sub_reports(self):
        """Create multiple mock sub-reports for testing."""
        return {
            "report1": MockSubReport("Report 1", "Content 1"),
            "report2": MockSubReport("Report 2", "Content 2"),
            "report3": MockSubReport("Report 3", "Content 3")
        }
    
    def test_reporter_initialization(self, reporter):
        """Test Reporter initialization with default sub-reports."""
        assert isinstance(reporter.sub_reports, dict)
        # Should have at least the default MutationEffectsReport
        assert len(reporter.sub_reports) > 0
        assert "Mutation Effects" in reporter.sub_reports
    
    def test_add_sub_report_valid(self, reporter):
        """Test adding a valid sub-report."""
        mock_report = MockSubReport("New Report", "New content")
        initial_count = len(reporter.sub_reports)
        
        reporter.add_sub_report(mock_report)
        
        assert len(reporter.sub_reports) == initial_count + 1
        assert "New Report" in reporter.sub_reports
        assert reporter.sub_reports["New Report"] == mock_report
    
    def test_add_sub_report_invalid_type(self, reporter):
        """Test adding invalid type raises TypeError."""
        with pytest.raises(TypeError, match="sub_report must be an instance of SubReport"):
            reporter.add_sub_report("not a sub-report")
    
    def test_remove_sub_report_exists(self, reporter):
        """Test removing an existing sub-report."""
        mock_report = MockSubReport("Removable Report", "Content")
        reporter.add_sub_report(mock_report)
        
        assert "Removable Report" in reporter.sub_reports
        reporter.remove_sub_report("Removable Report")
        assert "Removable Report" not in reporter.sub_reports
    
    def test_remove_sub_report_not_exists(self, reporter):
        """Test removing non-existent sub-report does nothing."""
        initial_count = len(reporter.sub_reports)
        reporter.remove_sub_report("Non-existent Report")
        assert len(reporter.sub_reports) == initial_count
    
    def test_list_sub_reports(self, reporter, mock_sub_reports):
        """Test listing all sub-report names."""
        # Add mock reports
        for report in mock_sub_reports.values():
            reporter.add_sub_report(report)
        
        report_names = reporter.list_sub_reports()
        
        assert isinstance(report_names, list)
        for name in mock_sub_reports.keys():
            # The names are the actual report names, not the keys
            report_name = mock_sub_reports[name].name
            assert report_name in report_names
    
    def test_generate_report_all_reports(self, reporter, sample_analyser, sample_mutations):
        """Test generating report with all sub-reports."""
        # Add some mock reports
        mock_report1 = MockSubReport("Mock Report 1", "Mock content 1")
        mock_report2 = MockSubReport("Mock Report 2", "Mock content 2")
        reporter.add_sub_report(mock_report1)
        reporter.add_sub_report(mock_report2)
        
        report = reporter.generate_report(sample_analyser, sample_mutations)
        
        # Should contain header
        assert "SNP ANALYSIS REPORT" in report
        # Should contain analyser info
        assert "Mitochondrial" in report
        assert str(sample_analyser.input_sequence) in report
        # Should contain mock report content
        assert "Mock content 1" in report
        assert "Mock content 2" in report
    
    def test_generate_report_specific_reports(self, reporter, sample_analyser, sample_mutations):
        """Test generating report with specific sub-reports."""
        mock_report1 = MockSubReport("Mock Report 1", "Mock content 1")
        mock_report2 = MockSubReport("Mock Report 2", "Mock content 2")
        reporter.add_sub_report(mock_report1)
        reporter.add_sub_report(mock_report2)
        
        # Generate report with only one sub-report
        report = reporter.generate_report(
            sample_analyser, 
            sample_mutations, 
            include_reports=["Mock Report 1"]
        )
        
        # Should contain the requested report
        assert "Mock content 1" in report
        # Should not contain the excluded report
        assert "Mock content 2" not in report
    
    def test_generate_report_unknown_reports(self, reporter, sample_analyser, sample_mutations):
        """Test generating report with unknown sub-report names raises ValueError."""
        with pytest.raises(ValueError, match="Unknown sub-reports"):
            reporter.generate_report(
                sample_analyser, 
                sample_mutations, 
                include_reports=["Non-existent Report"]
            )
    
    def test_generate_report_empty_include_list(self, reporter, sample_analyser, sample_mutations):
        """Test generating report with empty include list."""
        report = reporter.generate_report(
            sample_analyser, 
            sample_mutations, 
            include_reports=[]
        )
        
        # Should still have header and basic info
        assert "SNP ANALYSIS REPORT" in report
        assert "Mitochondrial" in report
        # Should not have any sub-report content
        lines = report.split('\n')
        content_lines = [line for line in lines if line.strip() and not line.startswith('=')]
        # Should be minimal content (just header info)
        assert len([line for line in content_lines if 'Mock' in line]) == 0
    
    def test_get_timestamp(self, reporter):
        """Test timestamp generation for report header."""
        timestamp = reporter._get_timestamp()
        
        # Should be in format YYYY-MM-DD HH:MM:SS
        assert len(timestamp) == 19  # Length of "2025-07-25 10:30:45"
        assert timestamp[4] == '-'  # Year separator
        assert timestamp[7] == '-'  # Month separator
        assert timestamp[10] == ' '  # Date/time separator
        assert timestamp[13] == ':'  # Hour separator
        assert timestamp[16] == ':'  # Minute separator
    
    def test_generate_report_includes_timestamp(self, reporter, sample_analyser, sample_mutations):
        """Test that generated report includes timestamp."""
        with patch.object(reporter, '_get_timestamp', return_value="2025-07-25 10:30:45"):
            report = reporter.generate_report(sample_analyser, sample_mutations)
            
            assert "Generated on: 2025-07-25 10:30:45" in report


class TestReporterIntegration:
    """Integration tests for Reporter with real sub-reports."""
    
    @pytest.fixture
    def reporter_with_mutation_effects(self):
        """Create a Reporter with the real MutationEffectsReport."""
        return Reporter()  # Already includes MutationEffectsReport by default
    
    @pytest.fixture
    def analyser_with_mutations(self):
        """Create an Analyser instance with applied mutations."""
        analyser = Analyser("ATGCGATCGTAG", dna_type="mitochondrial")
        mutations = [{"position": 3, "alt": "T"}]
        analyser.apply_mutations(mutations)
        return analyser, mutations
    
    def test_full_report_generation(self, reporter_with_mutation_effects, analyser_with_mutations):
        """Test generating a complete report with real components."""
        analyser, mutations = analyser_with_mutations
        
        report = reporter_with_mutation_effects.generate_report(analyser, mutations)
        
        # Should contain main header
        assert "SNP ANALYSIS REPORT" in report
        # Should contain analyser information
        assert "Mitochondrial" in report
        assert str(analyser.input_sequence) in report
        # Should contain mutation effects report
        assert "Mutation Effects" in report
        # Should have proper formatting throughout
        assert "=" * ReportFormatter.SECTION_WIDTH in report
    
    def test_multiple_mutations_report(self, reporter_with_mutation_effects):
        """Test report generation with multiple mutations."""
        analyser = Analyser("ATGCGATCGTAGCCG", dna_type="nuclear")
        mutations = [
            {"position": 3, "alt": "T"},
            {"position": 8, "alt": "A"},
            {"position": 12, "alt": "T"}
        ]
        analyser.apply_mutations(mutations)
        
        report = reporter_with_mutation_effects.generate_report(analyser, mutations)
        
        # Should contain information about all mutations
        assert "Number of Mutations: 3" in report
        # Should be formatted properly
        assert "Nuclear" in report
