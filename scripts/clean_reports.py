#!/usr/bin/env python3
"""
Clean reports database and generated reports directory.
Usage: python clean_reports.py [--dry-run]
"""

import os
import sys
import shutil
from pathlib import Path

def main():
    dry_run = "--dry-run" in sys.argv
    
    # Paths
    repo_root = Path(__file__).parent.parent
    db_path = repo_root / "backend" / "analyser_backend" / "report_jobs.sqlite3"
    reports_dir = repo_root / "backend" / "generated_reports"
    
    print("=" * 60)
    print("SNP Analyser: Clean Reports & Database")
    print("=" * 60)
    
    if dry_run:
        print("[DRY RUN] No changes will be made.\n")
    
    # 1. Database
    if db_path.exists():
        size_mb = db_path.stat().st_size / (1024 * 1024)
        print(f"\n📊 Database: {db_path}")
        print(f"   Size: {size_mb:.2f} MB")
        if not dry_run:
            db_path.unlink()
            print("   ✓ Deleted")
        else:
            print("   [DRY RUN] Would delete")
    else:
        print(f"\n📊 Database: {db_path}")
        print("   (does not exist)")
    
    # 2. Reports directory
    if reports_dir.exists():
        report_files = list(reports_dir.glob("*"))
        print(f"\n📁 Reports directory: {reports_dir}")
        print(f"   Files: {len(report_files)}")
        
        if report_files:
            for f in sorted(report_files)[:5]:  # Show first 5
                print(f"     - {f.name}")
            if len(report_files) > 5:
                print(f"     ... and {len(report_files) - 5} more")
        
        if not dry_run and report_files:
            shutil.rmtree(reports_dir)
            reports_dir.mkdir(exist_ok=True)
            print("   ✓ Deleted all reports")
        elif dry_run and report_files:
            print("   [DRY RUN] Would delete all reports")
    else:
        print(f"\n📁 Reports directory: {reports_dir}")
        print("   (does not exist)")
    
    print("\n" + "=" * 60)
    if dry_run:
        print("Run without --dry-run to actually delete:")
        print(f"  python {__file__}")
    else:
        print("✓ Cleanup complete")
    print("=" * 60)

if __name__ == "__main__":
    main()
