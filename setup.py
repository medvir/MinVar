#!/usr/bin/env python

import setuptools

setuptools.setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive',
                    'isort', 'pytest-runner'],
    install_requires=['setuptools_scm', 'biopython'],
    tests_require=['pytest', 'flake8'],
    name='MinVar',
    description='Minority variants in HIV',
    url='https://github.com/ozagordi/MinVar',
    author='Osvaldo Zagordi',
    author_email='firstname.lastname@gmail.com',
    packages=setuptools.find_packages('src'),
    scripts=['bin/minvar', 'src/scripts/blast2sam.py'],
    #package_dir={'minvar': 'src/minvar'},
    package_dir={'': 'src'},
    package_data={'minvar': ['db/*']},
    entry_points={
        'console_scripts': ['minvar = minvar.cli:main']
        }
)
