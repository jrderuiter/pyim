__author__ = 'Julian'

import sys

from setuptools import setup, find_packages

install_requires = [
    'numpy',
    'scipy',
    'pandas',
    'pysam',
    'natsort',
    'rpy2'
]

setup(
    name='pyim',
    version='0.3',
    url='',
    author='Julian de Ruiter',
    author_email='j.r.deruiter@icloud.com',
    description='Predicts transposon insertion sites from DNA-seq data.',
    license='BSD',
    packages=find_packages(),
    include_package_data=True,
    entry_points={'console_scripts': [
        'pyim-align = pyim.main.main_alignment:main',
        'pyim-annotate = pyim.main.main_annotation:main',
        'pyim-cis = pyim.main.main_cis:main'
    ]},
    extras_require={},
    zip_safe=True,
    classifiers=[],
    install_requires=install_requires
)
