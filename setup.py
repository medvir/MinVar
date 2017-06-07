#!/usr/bin/env python

import setuptools

setuptools.setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
    install_requires=['setuptools_scm'],
    name='MinVar',
      description='Minority variants in HIV',
      url='https://github.com/ozagordi/MinVar',
      author='Osvaldo Zagordi',
      author_email='firstname.lastname@gmail.com',
      packages = setuptools.find_packages(),
      scripts = ['bin/minvar', 'minvar/blast2sam.py'],
      data_files = [('minvar/db',
        ['minvar/db/HXB2_pol_gene.fasta',
         'minvar/db/HIV_consensus.fasta',
         'minvar/db/HIV_cons_db.nsq',
         'minvar/db/HIV_cons_db.nin',
         'minvar/db/HIV_cons_db.nhr',
         'minvar/db/template.tex',
         'minvar/db/masterComments_RTI.txt',
         'minvar/db/masterComments_PI.txt',
         'minvar/db/masterComments_INI.txt',
         'minvar/db/consensus_B.fna',
         'minvar/db/consensus_B.faa',
         'minvar/db/consensus_B.bed',
         'minvar/db/HXB2_pol_gene.fasta.fai',
         'minvar/db/consensus_B.fna.fai'])],
      entry_points={
          'console_scripts': [
              'minvar = minvar.__main__:main'
          ]
    }
)
