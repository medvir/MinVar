# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
- Bugfix: when no RAS mutation is found, `reportdrm` failed.
- Bugfix: when the frequency of a mutation is 0.5/0.5, `cns_max_freq` was longer than the reference.
- Problems with an assert when all reads have indels: removed the assert, `phase_mutations` has limitations.

## [v2.1.2] - 2017-12-19
- More information on the run included in the report.

## [v2.1.1] - 2017-12-19
- Fixed a bug showing up when sample consensus starts with a stop codon.

## [v2.1] - 2017-12-18

### Added
- Support for non-overlapping amplicons.

## [v2.0] - 2017-12-13

### Added
- Support for HCV.
- Report in pdf now styled as article with information on footer of each page.
- Contact information and a logo can be specified with a INI file in `~/.minvar`
- `minvar` can be called as a module or as command line program.

### Changed
- Improved generation of sample consensus.
- Widespread use of `shlex`, _i.e._ processes are not called via shell.
- Reorganization of the package structure:
  - code is now in `src/minvar`;
  - data are now in `src/minvar/db`;
  - data now include plenty of resources for HCV.
- Import system has been changed.
- Tests have been adapted (coverage is still very low).
- Code quality via flake8, pep8, isort and pep257.

## [v1.2b] - 2017-07-08
### Changed
- `sort` in `pandas` has been replaced by `sort_values`.
- Latex template needs to define `\highlight`.

## [v1.2a3] - 2017-06-13
### Added
- MinVar is now on bioconda.

## [v1.2alpha] - 2017-06-07
### Added
- `minvar -v` now returns version

## [v1.1] - 2017-03-08
### Added
- Switch `-r` added to turn on/off the calibration with GATK.
- Documentation.
- Ansible instructions.

## [v1.0.1] - 2016-09-20
Production version, tested on 454 and Miseq and used in the manuscript.
