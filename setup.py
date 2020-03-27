from setuptools import find_packages, setup, Command
import distutils.command.build
from distutils.dist import Distribution
# from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
# import xmlrunner
import unittest
import sys
import os
# import dotenv
import logging
import glob
from datetime import datetime
import input

# Build the Python package
setup(
    name='input',
    version=input.VERSION,
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'input=input.command_line:input_cli'
        ],
    },
    author='Franziska Lang',
    description='TODO',
    requires=[],
    # NOTE: always specify versions to ensure build reproducibility
    install_requires=['biopython==1.76'],
    setup_requires=[],
    classifiers=[
        'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 2.7'
      ]
)
