# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.4.0 - 2024-06-25
### Added
- Add vidarr labels to outputs (changes to medata only).
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) 

## 1.3.0 - 2023-11-28
### Added
- Introducing chromosome-size based coeffient for proportional allocation of job memory.

## 1.2.0 - 2022-06-10
### Changed
- Moved assembly-specific choices into the workflow.
- Updated README.md.
[GRD-635](https://jira.oicr.on.ca/browse/GRD-635)

## 1.1.3 - 2022-12-22
### Added
-Patching aggregate function so it handle empty inputs.

## 1.1.2 - 2022-12-23
### Added
- Introduced a parameter to set a minimum threshold for reads in input bam files.
- [GDI-2555](https://jira.oicr.on.ca/browse/GDI-2555)

## 1.1.1 - 2022-12-07
### Fixed
- Scattering by chromosome also brought in a possibility to encounter empty inputs in shards of extractMedipsCounts task. This is a fix to that.
- [GDI-2555](https://jira.oicr.on.ca/browse/GDI-2555)

## 1.1.0 - 2022-11-20
### Changed
- Scattering the task extractMedipsCounts by chromosome, updated default module for medips script. 

### Added
- Additional tasks for handling splitted outputs.
- [GDI-2615](https://jira.oicr.on.ca/browse/GDI-2615)

## 1.0.7 - 2022-06-07
### Changed
- Restoring reference argument for extract medips counts task.
- [GP-3386](https://jira.oicr.on.ca/browse/GP-3386)

## 1.0.6 - 2022-05-13
### Changed
- Making insert size metrics an optional output (to prevent fails when not enough data are available).
- [GP-2942](https://jira.oicr.on.ca/browse/GP-2942)

## 1.0.5 - 2022-04-19
### Fixed
- Small bug in the final shell script, to deal with grep causing pipefail.
- [GC-9245](https://jira.oicr.on.ca/browse/GC-9245)

### Added
- Making RT scripts more robust.
- [GP-2872](https://jira.oicr.on.ca/browse/GP-2872) 

## 1.0.4 - 2021-06-01
### Changed
- Migration to Vidarr.
