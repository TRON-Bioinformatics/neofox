from setuptools import find_packages, setup, Command
import neofox


# parses requirements from file
with open('requirements.txt') as f:
    required = f.read().splitlines()

# Build the Python package
setup(
    name='neofox',
    version=neofox.VERSION,
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'neofox=neofox.command_line:neofox_cli'
        ],
    },
    author=['Franziska Lang', 'Pablo Riesgo Ferreiro'],
    description='TODO',
    requires=[],
    # NOTE: always specify versions to ensure build reproducibility
    # NOTE2: sklearn==0.19.0 is a hidden dependency as it is required by Classifier.pickle
    install_requires=required,
    setup_requires=[],
    classifiers=[
        'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3'
      ]
)
