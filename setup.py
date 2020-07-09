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
    author=['Franziska Lang', 'Pablo Riesgo Ferreiro'],
    description='TODO',
    requires=[],
    # NOTE: always specify versions to ensure build reproducibility
    # NOTE2: sklearn==0.19.0 is a hidden dependency as it is required by Classifier.pickle
    install_requires=[
        'biopython==1.76',
        'mock',
        'pandas==0.24.2',
        'numpy==1.16.2',
        'scipy==1.4.1',
        'pickle-mixin',
        'scikit-learn==0.20.3',
        'logzero==1.5.0',
        'python-dotenv==0.12.0',
        'betterproto==1.2.5'
    ],
    setup_requires=[],
    classifiers=[
        'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3'
      ]
)
