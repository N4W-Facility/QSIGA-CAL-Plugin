# -*- coding: utf-8 -*-
"""
QSIGA-CAL Dependencies Auto-Installer

This module automatically installs missing Python packages when the plugin loads.

Author: QSIGA-CAL Development Team
"""

import subprocess
import sys
import os
from qgis.PyQt.QtWidgets import QMessageBox, QProgressDialog
from qgis.PyQt.QtCore import Qt
from qgis.core import QgsMessageLog, Qgis


# List of required packages
# Note: numpy 1.26.4 and scipy 1.13.0 are specific versions required for compatibility
REQUIRED_PACKAGES = {
    'numpy': 'numpy==1.26.4',
    'pandas': 'pandas>=1.1.0',
    'scipy': 'scipy==1.13.0',
    'spotpy': 'spotpy>=1.5.14',
    'matplotlib': 'matplotlib>=3.3.0',
    'openpyxl': 'openpyxl>=3.0.0',
    'xlrd': 'xlrd>=1.2.0',
    'rasterio': 'rasterio>=1.2.0',
    'numba': 'numba>=0.55.0',  # Required by pysheds, needs numpy<2.0
    'pysheds': 'pysheds>=0.3.0',
    'geopandas': 'geopandas>=0.9.0',
    'shapely': 'shapely>=1.7.0',
    'pyproj': 'pyproj>=3.0.0',
    'fiona': 'fiona>=1.8.0',
    'requests': 'requests>=2.25.0',
    'statsmodels': 'statsmodels>=0.12.0'
}


def get_python_executable():
    """
    Get the correct Python executable for QGIS.

    Returns:
        str: Path to Python executable
    """
    # On Windows, sys.executable might point to qgis-bin.exe
    # We need to find python.exe or python3.exe

    if sys.platform == 'win32':
        # Try to find Python in QGIS installation
        qgis_exe = sys.executable
        qgis_dir = os.path.dirname(qgis_exe)

        # Common paths in QGIS installation
        possible_paths = [
            os.path.join(qgis_dir, 'python.exe'),
            os.path.join(qgis_dir, 'python3.exe'),
            os.path.join(qgis_dir, '..', 'bin', 'python.exe'),
            os.path.join(qgis_dir, '..', 'bin', 'python3.exe'),
            # OSGeo4W paths
            'C:\\OSGeo4W\\bin\\python3.exe',
            'C:\\OSGeo4W64\\bin\\python3.exe',
        ]

        # Also check Program Files
        program_files = os.environ.get('PROGRAMFILES', 'C:\\Program Files')
        for version in ['3.44.4', '3.34', '3.32', '3.30', '3.28']:
            possible_paths.extend([
                os.path.join(program_files, f'QGIS {version}', 'bin', 'python3.exe'),
                os.path.join(program_files, f'QGIS {version}', 'bin', 'python.exe'),
            ])

        # Find the first existing Python executable
        for path in possible_paths:
            if os.path.exists(path):
                return path

        # Fallback: try using python from PATH
        return 'python'
    else:
        # On Linux/Mac, sys.executable should work fine
        return sys.executable


def check_missing_packages():
    """
    Check which required packages are missing or have incompatible versions.

    Returns:
        list: List of missing/incompatible package names
    """
    missing = []

    for package_name in REQUIRED_PACKAGES.keys():
        try:
            module = __import__(package_name)

            # Special check for numpy version (must be <2.0 for numba compatibility)
            if package_name == 'numpy':
                if hasattr(module, '__version__'):
                    version = module.__version__
                    major_version = int(version.split('.')[0])
                    if major_version >= 2:
                        # numpy 2.x is incompatible with numba
                        QgsMessageLog.logMessage(
                            f"Incompatible numpy version detected: {version} (need <2.0 for numba)",
                            "QSIGA-CAL",
                            Qgis.Warning
                        )
                        missing.append(package_name)
        except ImportError:
            missing.append(package_name)

    return missing


def install_packages(packages, parent=None):
    """
    Install missing packages using pip.

    Args:
        packages: List of package names to install
        parent: Parent widget for dialogs

    Returns:
        tuple: (success: bool, failed_packages: list)
    """
    if not packages:
        return True, []

    failed = []

    # Create progress dialog
    progress = QProgressDialog(
        "Installing missing dependencies...",
        "Cancel",
        0,
        len(packages),
        parent
    )
    progress.setWindowTitle("QSIGA-CAL Setup")
    progress.setWindowModality(Qt.WindowModal)
    progress.show()

    for i, package_name in enumerate(packages):
        if progress.wasCanceled():
            failed.extend(packages[i:])
            break

        progress.setLabelText(f"Installing {package_name}...")
        progress.setValue(i)

        package_spec = REQUIRED_PACKAGES[package_name]

        try:
            QgsMessageLog.logMessage(
                f"Installing {package_spec}...",
                "QSIGA-CAL",
                Qgis.Info
            )

            # Install package using correct Python executable
            python_exe = get_python_executable()

            QgsMessageLog.logMessage(
                f"Using Python: {python_exe}",
                "QSIGA-CAL",
                Qgis.Info
            )

            # Don't capture output to avoid issues with sys.stderr/stdout being None
            # The output will go to the QGIS log automatically

            # Force reinstall for numpy to ensure correct version
            extra_args = []
            if package_name == 'numpy':
                extra_args = ["--force-reinstall", "--no-deps"]  # Reinstall numpy without dependencies

            result = subprocess.run(
                [python_exe, "-m", "pip", "install", package_spec, "--quiet"] + extra_args,
                timeout=300,  # 5 minutes timeout
                stderr=subprocess.DEVNULL,  # Suppress stderr to avoid numpy issues
                stdout=subprocess.DEVNULL   # Suppress stdout
            )

            if result.returncode != 0:
                QgsMessageLog.logMessage(
                    f"Failed to install {package_name} (exit code: {result.returncode})",
                    "QSIGA-CAL",
                    Qgis.Warning
                )
                failed.append(package_name)
            else:
                QgsMessageLog.logMessage(
                    f"Successfully installed {package_name}",
                    "QSIGA-CAL",
                    Qgis.Success
                )

        except subprocess.TimeoutExpired:
            QgsMessageLog.logMessage(
                f"Timeout installing {package_name}",
                "QSIGA-CAL",
                Qgis.Warning
            )
            failed.append(package_name)
        except Exception as e:
            QgsMessageLog.logMessage(
                f"Error installing {package_name}: {str(e)}",
                "QSIGA-CAL",
                Qgis.Warning
            )
            failed.append(package_name)

    progress.setValue(len(packages))
    progress.close()

    return len(failed) == 0, failed


def ensure_dependencies(parent=None, auto_install=True):
    """
    Check and optionally install missing dependencies.

    Args:
        parent: Parent widget for dialogs
        auto_install: If True, automatically install missing packages

    Returns:
        bool: True if all dependencies are available
    """
    missing = check_missing_packages()

    if not missing:
        QgsMessageLog.logMessage(
            "All dependencies are installed",
            "QSIGA-CAL",
            Qgis.Info
        )
        return True

    QgsMessageLog.logMessage(
        f"Missing packages: {', '.join(missing)}",
        "QSIGA-CAL",
        Qgis.Warning
    )

    if not auto_install:
        # Show error message
        QMessageBox.critical(
            parent,
            "Missing Dependencies",
            f"The following required packages are missing:\n\n"
            f"{', '.join(missing)}\n\n"
            f"Please install them manually using:\n"
            f"python -m pip install {' '.join(REQUIRED_PACKAGES[p] for p in missing)}"
        )
        return False

    # Ask user if they want to install
    reply = QMessageBox.question(
        parent,
        "Install Missing Dependencies",
        f"QSIGA-CAL requires the following packages:\n\n"
        f"{chr(10).join('• ' + p for p in missing)}\n\n"
        f"Do you want to install them automatically?\n\n"
        f"(This may take a few minutes)",
        QMessageBox.Yes | QMessageBox.No,
        QMessageBox.Yes
    )

    if reply == QMessageBox.No:
        QMessageBox.warning(
            parent,
            "Dependencies Required",
            "QSIGA-CAL requires these packages to function properly.\n\n"
            "You can install them manually using:\n"
            f"python -m pip install {' '.join(REQUIRED_PACKAGES[p] for p in missing)}"
        )
        return False

    # Install packages
    success, failed = install_packages(missing, parent)

    if success:
        QMessageBox.information(
            parent,
            "Installation Complete",
            "All dependencies have been installed successfully!\n\n"
            "Please restart QGIS for the changes to take effect."
        )
        return True
    else:
        QMessageBox.warning(
            parent,
            "Installation Incomplete",
            f"Failed to install some packages:\n\n"
            f"{', '.join(failed)}\n\n"
            f"Please install them manually using:\n"
            f"python -m pip install {' '.join(REQUIRED_PACKAGES[p] for p in failed)}"
        )
        return False
