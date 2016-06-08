import sys

import setuptools
import versioneer

INSTALL_REQUIRES = ['future', 'numpy', 'scipy', 'pandas', 'pysam',
                    'rpy2', 'scikit-bio', 'toolz', 'tqdm', 'intervaltree']

EXTRAS_REQUIRE = {
    'dev': ['sphinx', 'pytest', 'pytest-mock',
            'pytest-datafiles', 'pytest-cov',
            'pytest-helpers-namespace']
}


# Check setuptools version, as recommended by:
# https://hynek.me/articles/conditional-python-dependencies/.
if int(setuptools.__version__.split('.', 1)[0]) < 18:
    assert 'bdist_wheel' not in sys.argv

    # Add pathlib for Pythons before 3.4.
    if sys.version_info[0:2] < (3, 4):
        INSTALL_REQUIRES.append('pathlib2')
else:
    EXTRAS_REQUIRE[":python_version<'3.4'"] = ['pathlib2']


setuptools.setup(
    name='pyim',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url='https://bitbucket.org/jrderuiter/pyim',
    author='Julian de Ruiter',
    author_email='julianderuiter@gmail.com',
    description='Predicts transposon insertion sites from DNA-seq data.',
    license='BSD',
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
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
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    zip_safe=True,
    classifiers=[]
)
