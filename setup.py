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
    version='0.4.0',
    url='',
    author='Julian de Ruiter',
    author_email='j.r.deruiter@icloud.com',
    description='Predicts transposon insertion sites from DNA-seq data.',
    license='BSD',
    packages=find_packages(),
    include_package_data=True,
    entry_points={'console_scripts': [
        'pyim-align = pyim.tools.alignment.main:main',
        'pyim-annotate = pyim.tools.annotation.main:main',
        'pyim-cis = pyim.tools.cis.main:main'
    ]},
    extras_require={},
    zip_safe=True,
    classifiers=[],
    install_requires=install_requires
)
