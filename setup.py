#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['pyfaidx', 'intervaltree', 'tqdm', 'toolz', 'rpy2', 'numpy',
                'pandas', 'pysam']

test_requirements = ['pytest', 'pytest-cov', 'pytest-mock',
                     'pytest-helpers-namespace', 'python-coveralls']

setup(
    name='pyim',
    version='0.2.0',
    description="Tools for analyzing insertional mutagenesis data.",
    long_description=readme + '\n\n' + history,
    author="Julian de Ruiter",
    author_email='julianderuiter@gmail.com',
    url='https://github.com/jrderuiter/pyim',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='pyim',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    extras_require={
        'test': test_requirements
    },
    entry_points={'console_scripts': [
        'pyim-align = pyim.main.pyim_align:main',
        'pyim-demultiplex = pyim.main.pyim_demultiplex:main',
        'pyim-merge = pyim.main.pyim_merge:main',
        'pyim-cis = pyim.main.pyim_cis:main',
        'pyim-annotate = pyim.main.pyim_annotate:main',
        'pyim-bed = pyim.main.pyim_bed:main',
        'pyim-split = pyim.main.pyim_split:main'
    ]})
