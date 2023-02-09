#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#
from setuptools import find_packages, setup, Command
import neofox


# parses requirements from file
with open("requirements.txt", encoding="utf-8") as f:
    required = f.read().splitlines()

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Build the Python package
setup(
    name="neofox",
    version=neofox.VERSION,
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "neofox=neofox.command_line:neofox_cli",
            "neofox-epitope=neofox.command_line:neofox_epitope_cli",
            "neofox-configure=neofox.command_line:neofox_configure",
        ],
    },
    author_email="franziska.lang@tron-mainz.de",
    author="TRON - Translational Oncology at the University Medical Center of the Johannes Gutenberg University Mainz "
    "- Computational Medicine group",
    description="Annotation of mutated peptide sequences (mps) with published or novel potential neo-epitope "
    "descriptors",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tron-bioinformatics/neofox",
    requires=[],
    # NOTE: always specify versions to ensure build reproducibility
    # NOTE2: sklearn==0.19.0 is a hidden dependency as it is required by Classifier.pickle
    install_requires=required,
    setup_requires=[],
    classifiers=[
        "Development Status :: 5 - Production/Stable",  # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Unix"
    ],
    python_requires='>=3.6,<3.9',
    license='GPLv3',
)
