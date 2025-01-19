# Changelog

## [1.5.2] - 2024-12-07

### Removed

- Obsolete `base-threads` dependency.

## [1.5.1] - 2024-11-23

### Changed

- Reformatted everything with ocamlformat.
- Minor fixes.

## [1.5.0] - 2019-11-22

### Changed

- Switched to OPAM file generation via `dune-project`.
- Updates for new Core API `v0.13.0`.

## [1.4.1] - 2018-10-24

### Changed

- Updated to OPAM 2.0.

## [1.4.0] - 2018-08-19

### Fixed

- Lacaml 11.0.0 changes.

### Changed

- Switched to dune and dune-release.

## [1.3.1] - 2017-08-05

### Fixed

- Incorrect threads dependency in OPAM package.

## [1.3.0] - 2017-07-30

### Changed

- Switched to jbuilder and topkg.

## Changes Before Version 1.3.0

```text
2017-05-23:  Fixes for Lacaml 9.2.3.

2017-01-12:  Switched to GPL 2 for the application.

2015-07-03:  Fixes wrt. new Lacaml API.

2014-10-03:  Using module aliases for improved compilation and linking time,
             and executable size.

2013-01-19:  Fixed build problem with newest Lacaml release.

2012-07-20:  Downgraded findlib version constraint to support the Debian
             testing branch.

2012-07-15:  New major release version 1.0.0:

               * Upgraded to OCaml 4.00
               * Upgraded to Core library 108.00.01
               * Switched to Oasis for packaging
               * Switched to OCamlBuild for the build process
               * Switched to GPL version 3 for application and
                 LGPL-2.1 with linking exception for library
               * Rewrote README in Markdown
               * Added stricter compilation flags

2012-02-01:  Updated to Lacaml 6.0.0.

2011-09-16:  Updated to Core 107.01.

2010-08-28:  Added support for hyper variables in inputs.

2010-08-14:  Fixed command-line help.

2010-07-11:  Added some support for online learning.

2009-11-23:  Initial public release.
```
