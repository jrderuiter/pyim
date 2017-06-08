#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    README = readme_file.read()

with open('HISTORY.rst') as history_file:
    HISTORY = history_file.read()

REQUIREMENTS = [
    'pyfaidx>=0.4.8.1', 'intervaltree>=2.1', 'tqdm>=4.7', 'toolz>=0.8',
    'rpy2>=2.8.2', 'numpy', 'pandas>=0.18', 'pysam>=0.9', 'natsort',
    'cutadapt >=1.8'
]

EXTRAS_REQUIRE = {
    'dev': [
        'pytest', 'pytest-cov', 'pytest-mock', 'pytest-helpers-namespace',
        'python-coveralls', 'sphinx', 'sphinx-autobuild', 'sphinx_rtd_theme',
        'bumpversion'
    ]
}

setup(
    name='pyim',
    version='0.2.1',
    description=('Tool for identifying transposon insertions '
                 'from targeted DNA-sequencing data.'),
    long_description=README + '\n\n' + HISTORY,
    author='Julian de Ruiter',
    author_email='julianderuiter@gmail.com',
    url='https://github.com/jrderuiter/pyim',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=REQUIREMENTS,
    license='MIT license',
    zip_safe=False,
    keywords='pyim',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    extras_require=EXTRAS_REQUIRE,
    entry_points={
        'console_scripts': [
            'pyim-align = pyim.main.pyim_align:main',
            'pyim-demultiplex = pyim.main.pyim_demultiplex:main',
            'pyim-merge = pyim.main.pyim_merge:main',
            'pyim-cis = pyim.main.pyim_cis:main',
            'pyim-annotate = pyim.main.pyim_annotate:main',
            'pyim-bed = pyim.main.pyim_bed:main',
            'pyim-split = pyim.main.pyim_split:main'
        ]
    })
