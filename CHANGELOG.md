# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Support for HCV.
- `minvar` can be called as a module or as command line program.

### Changed
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
