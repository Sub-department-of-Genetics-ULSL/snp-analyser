"""
Base classes for the reporting system.

This module contains the abstract base classes and utilities that are shared
across the reporting system to avoid circular imports.
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional


class ReportFormatter:
    """Utility class for consistent report formatting."""
    
    # Standard formatting constants
    SECTION_WIDTH = 80
    SUBSECTION_WIDTH = 60
    DETAIL_WIDTH = 40
    
    # Standard characters
    SECTION_CHAR = "="
    SUBSECTION_CHAR = "-"
    DETAIL_CHAR = "."
    
    @classmethod
    def format_section_header(cls, title: str) -> str:
        """Format a main section header with consistent boundaries."""
        padded_title = f" {title.upper()} "
        padding = (cls.SECTION_WIDTH - len(padded_title)) // 2
        return f"{cls.SECTION_CHAR * padding}{padded_title}{cls.SECTION_CHAR * padding}"
    
    @classmethod
    def format_subsection_header(cls, title: str) -> str:
        """Format a subsection header with consistent boundaries."""
        padded_title = f" {title} "
        padding = (cls.SUBSECTION_WIDTH - len(padded_title)) // 2
        # Ensure total width is exactly SUBSECTION_WIDTH
        remaining = cls.SUBSECTION_WIDTH - len(padded_title) - (2 * padding)
        return f"{cls.SUBSECTION_CHAR * padding}{padded_title}{cls.SUBSECTION_CHAR * (padding + remaining)}"
    
    @classmethod
    def format_section_separator(cls) -> str:
        """Create a section separator line."""
        return cls.SECTION_CHAR * cls.SECTION_WIDTH
    
    @classmethod
    def format_subsection_separator(cls) -> str:
        """Create a subsection separator line."""
        return cls.SUBSECTION_CHAR * cls.SUBSECTION_WIDTH
    
    @classmethod
    def format_detail_separator(cls) -> str:
        """Create a detail separator line."""
        return cls.DETAIL_CHAR * cls.DETAIL_WIDTH
    
    @classmethod
    def format_key_value(cls, key: str, value: str, indent: int = 0) -> str:
        """Format key-value pairs consistently."""
        spaces = " " * indent
        return f"{spaces}{key}: {value}"
    
    @classmethod
    def format_list_item(cls, item: str, bullet: str = "•", indent: int = 2) -> str:
        """Format list items consistently."""
        spaces = " " * indent
        return f"{spaces}{bullet} {item}"
    
    @classmethod
    def wrap_content_block(cls, content: str, title: Optional[str] = None) -> str:
        """Wrap content in a consistent block format."""
        lines = []
        
        if title:
            lines.append(cls.format_subsection_header(title))
        else:
            lines.append(cls.format_subsection_separator())
        
        lines.append("")
        lines.append(content)
        lines.append("")
        lines.append(cls.format_subsection_separator())
        
        return "\n".join(lines)


class SubReport(ABC):
    """Abstract base class for sub-reports."""
    
    @abstractmethod
    def generate_content(self, analyser, mutations: List[Dict]) -> str:
        """Generate the main content of the sub-report.
        
        Args:
            analyser: Analyser instance with sequence data
            mutations: List of mutations that were applied
            
        Returns:
            String containing the sub-report content (without boundaries)
        """
        pass
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Return the name of this sub-report."""
        pass
    
    def generate(self, analyser, mutations: List[Dict]) -> str:
        """Generate a formatted sub-report with consistent boundaries.
        
        Args:
            analyser: Analyser instance with sequence data
            mutations: List of mutations that were applied
            
        Returns:
            String containing the formatted sub-report
        """
        content = self.generate_content(analyser, mutations)
        return ReportFormatter.wrap_content_block(content, self.name)
