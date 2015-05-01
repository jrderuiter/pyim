__author__ = 'Julian'

import sys

from setuptools import setup, find_packages

install_requires = ['future', 'numpy', 'scipy', 'pandas', 'pysam',
                    'natsort', 'rpy2', 'scikit-bio', 'tkgeno']

if sys.version_info[0] == 2:
    install_requires += ['pathlib', 'enum']

setup(
    name='pyim',
    version='0.4.2',
    url='',
    author='Julian de Ruiter',
    author_email='j.r.deruiter@icloud.com',
    description='Predicts transposon insertion sites from DNA-seq data.',
    license='BSD',
    packages=find_packages(),
    include_package_data=True,
    entry_points={'console_scripts': [
        'pyim-align = pyim.main.align:main',
        'pyim-merge = pyim.main.merge:main',
        'pyim-annotate = pyim.main.annotate:main',
        'pyim-cis = pyim.tools.cis.main:main'
    ]},
    extras_require={'test': 'pytest'},
    zip_safe=True,
    classifiers=[],
    install_requires=install_requires
)
