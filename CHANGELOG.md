# Changelog

All notable changes to this project will be documented in this file.

## [1.1.0] - 2026-01-05

### Added

- **CRYSTAL23 support** (`.d12`) with linting, hover, completion, and semantic tokens
- **Formatting** for INCAR, POSCAR/CONTCAR, and KPOINTS
- **Custom command** `vasp.previewStructure` for AST inspection
- **Filetype overrides** via LSP settings and initialization options
- **Release automation** on tag push (build, test, npm pack + GitHub Release)

### Changed

- **Document cache** eviction upgraded to LRU
- **KPOINTS semantic tokens** now respect explicit vs automatic sections
- **POTCAR validation** falls back to CONTCAR when POSCAR is missing
- Documentation refreshed across user and contributor guides

### Fixed

- Diagnostic consistency across POSCAR parsing
- Stress‑test thresholds for large POSCAR files in local environments

## [1.0.0] - 2025-12-25

### Added

- **Modular Project Structure**: Feature-based organization (`src/features/`)
- **Semantic Tokens**: Rich highlighting for VASP INCAR files
- **Folding Ranges**: POSCAR lattice/coordinate folding + comments
- **Code Actions**: Quick fixes for INCAR tag typos
- **Document Symbols**: Outline view for POSCAR/CONTCAR + INCAR
- **Hover Documentation**: VASP Wiki tag docs
- **Auto-completion**: KPOINTS snippets + INCAR tag suggestions
- **Cross-file Linting**: POTCAR order validation vs POSCAR
- **CI Pipeline**: Lint, build, and test in GitHub Actions

### Fixed

- File URI → path conversion for cross-file validation
- ESLint and Prettier warnings
- Duplicate key errors in `package.json`
- Missing build script in `package.json`

### Changed

- Refactored entry point into a class-based `LspServer`
- Documentation updates with Neovim setup guide

---
*Generated for the VASP Language Server Launch*
