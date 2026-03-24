#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QSIGA-CAL Plugin Packager for QGIS Installation

Creates a single ZIP file of the complete QSIGA_CAL plugin that can be installed
directly from QGIS using: Plugins → Install from ZIP

Usage:
    python package_qsiga_cal.py

Output:
    QSIGA_CAL_v0.1_YYYYMMDD.zip (ready for QGIS installation)

Author: QSIGA-CAL Development Team
"""

import os
import zipfile
from datetime import datetime


def should_exclude(file_path, exclude_patterns):
    """Check if a file should be excluded from the package."""
    file_path = file_path.replace('\\', '/')
    for pattern in exclude_patterns:
        if pattern in file_path:
            return True
    return False


def create_qsiga_cal_zip(source_dir, output_dir):
    """
    Create a ZIP file of QSIGA_CAL plugin for QGIS installation.

    Args:
        source_dir: Path to qsiga_cal directory
        output_dir: Where to save the ZIP file
    """
    print("="*80)
    print("QSIGA-CAL PLUGIN PACKAGER FOR QGIS")
    print("="*80)
    print()

    # Get version from metadata
    version = "0.1"
    metadata_path = os.path.join(source_dir, 'metadata.txt')
    if os.path.exists(metadata_path):
        with open(metadata_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('version='):
                    version = line.split('=')[1].strip()
                    break

    # Create output filename
    timestamp = datetime.now().strftime("%Y%m%d")
    output_filename = f"QSIGA_CAL_v{version}_{timestamp}.zip"
    output_path = os.path.join(output_dir, output_filename)

    print(f"Plugin version: {version}")
    print(f"Output file: {output_path}")
    print()

    # Files and directories to exclude
    exclude_patterns = [
        '__pycache__',
        '.pyc',
        '.git/',
        '.gitignore',
        '.vscode/',
        '.idea/',
        '.claude/',
        'ejemplo_funciona',  # Old working example
        'Step_03/Data/',  # Test data in Step_03 only (exclude calibration test files)
        'temp/',
        'temp_unified_plugin/',
        '.DS_Store',
        'Thumbs.db',
        '.pytest_cache',
        'pb_tool.cfg',
        'pylintrc',
        'Makefile',
        'plugin_upload.py',
        'package_plugins.py',
        'package_for_qgis.py',
        'package_qsiga_cal.py',  # Don't include this packaging script
        'fix_subplugins.py',  # Don't include fix script
        'create_unified_plugin.py',
        'INSTRUCCIONES_EMPAQUETADO.md',
        'qsiga_plugins_for_qgis/',
        'earthdata_token.txt',
        '/test/',
        'ejecutar_siga_simple.py',
        'create_timeseries_siga.py',
        'fix_timeseries_qcapt.py',
        'process_lai.py',
        'metadata.txt.bak'  # Exclude backup metadata files (keep them renamed)
    ]

    print("Creating package...")
    print("-" * 80)
    print()

    files_added = 0
    files_excluded = 0

    with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(source_dir):
            # Remove excluded directories
            dirs[:] = [d for d in dirs if not should_exclude(
                os.path.join(root, d), exclude_patterns)]

            for file in files:
                file_path = os.path.join(root, file)

                # Check if file should be excluded
                if should_exclude(file_path, exclude_patterns):
                    files_excluded += 1
                    continue

                # Create archive name (relative path from parent of source_dir)
                # IMPORTANT: Must preserve QSIGA_CAL folder name as root in ZIP
                rel_path = os.path.relpath(file_path, os.path.dirname(source_dir))
                arcname = rel_path.replace('qsiga_cal', 'QSIGA_CAL')  # Ensure proper casing

                # Add file to zip
                zipf.write(file_path, arcname)
                files_added += 1

                # Show progress
                if files_added % 50 == 0:
                    print(f"  Added {files_added} files...")

    size_mb = os.path.getsize(output_path) / (1024 * 1024)

    print()
    print("="*80)
    print("PACKAGE CREATED SUCCESSFULLY!")
    print("="*80)
    print()
    print(f"Package location: {output_path}")
    print(f"Files added: {files_added}")
    print(f"Files excluded: {files_excluded}")
    print(f"Package size: {size_mb:.2f} MB")
    print()
    print("="*80)
    print("INSTALLATION INSTRUCTIONS")
    print("="*80)
    print()
    print("STEP 1: Install Dependencies")
    print("-" * 80)
    print("Before installing the plugin, install required Python packages:")
    print()
    print("Windows (run as Administrator):")
    print("  1. Open OSGeo4W Shell")
    print("  2. Run: python -m pip install numpy pandas scipy spotpy matplotlib openpyxl xlrd")
    print()
    print("Or use the install_dependencies.bat file from the Plugins folder.")
    print()
    print("STEP 2: Install Plugin in QGIS")
    print("-" * 80)
    print("1. Open QGIS")
    print("2. Go to: Plugins → Manage and Install Plugins")
    print("3. Click on: 'Install from ZIP'")
    print(f"4. Select: {output_filename}")
    print("5. Click: 'Install Plugin'")
    print("6. Wait for confirmation message")
    print()
    print("STEP 3: Activate Plugin")
    print("-" * 80)
    print("1. In QGIS: Plugins → Manage and Install Plugins")
    print("2. Go to 'Installed' tab")
    print("3. Check the box for 'QSIGA-CAL'")
    print("4. Close the dialog")
    print()
    print("The plugin menu 'QSIGA-CAL' will appear in QGIS with all 5 steps!")
    print()

    return output_path


if __name__ == "__main__":
    # Get script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Source is the current directory (qsiga_cal)
    source_dir = script_dir

    # Output to parent directory
    output_dir = os.path.dirname(script_dir)

    print()
    print("This script creates a single ZIP package of the complete QSIGA-CAL plugin")
    print("that can be installed directly from QGIS.")
    print()
    print(f"Source directory: {source_dir}")
    print(f"Output directory: {output_dir}")
    print()

    response = input("Continue? (y/n): ").strip().lower()

    if response == 'y' or response == 'yes':
        print()
        create_qsiga_cal_zip(source_dir, output_dir)
        print()
    else:
        print("Operation cancelled.")
