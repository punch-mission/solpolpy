# Changelog

[Available through GitHub](https://github.com/punch-mission/solpolpy/releases)

## Version 0.4.1: July 22, 2025

- Add RGB image visualization and sped up resolution in https://github.com/punch-mission/solpolpy/pull/185
- Allows selecting a non-default HDU in https://github.com/punch-mission/solpolpy/pull/189
- Adds docs button to view source in https://github.com/punch-mission/solpolpy/pull/190

## Version 0.4.0: Nov 22, 2024

- Fix IMAX effect to use offset by @jmbhughes in #175
- Fix nfi imax distortion by @s0larish in #176
- This release splits the "mzp" system into "mzpsolar" and "mzpinstru" to represent the solar and instrument referenced systems separately.

## Version 0.3.4: Nov 13, 2024

- add offset angle to resolve function by @jmbhughes in #173

## Version 0.3.3: Nov 2, 2024

- avoid long versions by @jmbhughes in #170
- fix error in npol_to_mzp transformation by @s0larish in #158

## Version 0.3.2: Nov 1, 2024

- check mzp to npol conversion with tests by @jmbhughes in #151
- Create CITATION.cff by @jmbhughes in #153
- add 3.13 to ci by @jmbhughes in #164
- Create CI_fixed.yaml by @jmbhughes in #167
- Switch to ReadTheDocs and versioning to automatic by @jmbhughes in #169

## Version 0.3.1: Aug 6, 2024

- fixes mzp to npol conversion by @jmbhughes in #149

## Version 0.3.0: Aug 6, 2024

- Use units everywhere! by @jmbhughes in #144

## Version 0.2.0: Jul 19, 2024

- Fixes doc script by @jmbhughes in #127
- [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in #128
- Add support for arbitrary polarization angles specification and npol system by @jmbhughes in #138

## Version 0.1.3: Jun 3, 2024

- added short description on IMAX effect by @s0larish in #112
- adds development guide note by @jmbhughes in #113
- bump python version to 3.10 by @jmbhughes in #125

## Version 0.1.2: Mar 8, 2024

### New features

This release allows for calculating the IMAX effect when doing conversions.

### What's Changed

- Update LICENSE.txt by @jmbhughes in #97
- Update citation by @jmbhughes in #102
- IMAX effect integration by @lowderchris in #98

## Version 0.1.1: Jan 12, 2024

This release is primarily to make sure masks are combined properly when doing conversions.

- Update cite.rst by @jmbhughes in #63
- Update citation by @jmbhughes in #64
- adds pre-commit by @jmbhughes in #65
- Adds pre-commit by @jmbhughes in #66
- [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in #75
- Update requirements.txt by @jmbhughes in #76
- Adds pre-commit, updates requirements by @jmbhughes in #77
- Adds test for angles in load_data, reformats some polarizers by @s0larish in #78
- Fix citation by @jmbhughes in #79
- moves citation, adds angle name test by @github-actions in #80
- added test for BtBr to MZP directly and indirectly by @jmbhughes in #81
- [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in #83
- [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in #85
- Weekly merge to develop by @github-actions in #84
- [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in #86
- Added image mask to carry around conversions #55 by @s0larish in #82
- Update version by @jmbhughes in #88
- fixes mask propagation, releases new version by @github-actions in #87

## Version 0.1.0: Dec 7, 2023

This release greatly improves the documentation. It also improves the resolve function while adding a new plot_collection function.

- Update README.md by @jmbhughes in #46
- Update readme by @jmbhughes in #47
- add doi badge by @jmbhughes in #49
- Adds DOI badge by @github-actions in #52
- v0.1.0 by @jmbhughes in #53
- Releases v0.1.0 by @jmbhughes in #56
- Fix docs deploy by @jmbhughes in #58
- test pandoc github action deploy by @jmbhughes in #59
- Add pandoc by @jmbhughes in #60
- add-pandoc by @jmbhughes in #61
- add pandoc setup by @jmbhughes in #62

## Version 0.0.1: Nov 6, 2023

- Universal pipeline by @jmbhughes in #13
- Added basic alpha support by @jmbhughes in #15
- Implemented astropy units support by @jmbhughes in #16
- added unit tests for sanitization by @jmbhughes in #17
- copied Bryce's tests and started incorporating by @jmbhughes in #19
- FITS functionality added by @jmbhughes in #21
- move instrument specific things to separate module by @jmbhughes in #22
- Transition to NDCube, Completed all equations by @jmbhughes in #23
- Create dependabot.yml by @jmbhughes in #36
- Create CODE_OF_CONDUCT.md by @jmbhughes in #38
- Update README.md by @jmbhughes in #39
- Create python-publish.yaml by @jmbhughes in #43
- pre-release updates by @s0larish in #42
- Create first release of solpolpy by @jmbhughes in #37
