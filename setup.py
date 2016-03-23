import sys

from setuptools import setup, find_packages

from version import get_git_version


install_requires = ['future', 'numpy', 'scipy', 'pandas', 'pysam',
                    'rpy2', 'scikit-bio', 'toolz', 'tqdm', 'intervaltree']

if not sys.version_info >= (3, ):
    install_requires += ['pathlib']

setup(
    name='pyim',
    version=get_git_version(),
    url='https://bitbucket.org/jrderuiter/pyim',
    author='Julian de Ruiter',
    author_email='julianderuiter@gmail.com',
    description='Predicts transposon insertion sites from DNA-seq data.',
    license='BSD',
    packages=find_packages(),
    include_package_data=True,
    entry_points={'console_scripts': [
        'pyim-align = pyim.main.align:main',
        'pyim-merge = pyim.main.merge:main',
        'pyim-annotate = pyim.main.annotate:main',
        'pyim-cis = pyim.main.cis:main',
        'pyim-plot = pyim.main.plot:main',
        'pyim-gff = pyim.main.gff:main',
        'pyim-split = pyim.main.split:main'
    ]},
    extras_require={'test': 'pytest'},
    zip_safe=True,
    classifiers=[],
    install_requires=install_requires
)
