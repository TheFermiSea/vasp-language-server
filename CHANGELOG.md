# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] - Launch Release - 2025-12-25

### Added

- **Modular Project Structure**: Refactored codebase to use a feature-based directory structure (`src/features/`).
- **Semantic Tokens**: Rich syntax highlighting for VASP INCAR files.
- **Folding Ranges**: Support for logical block folding in POSCAR (Lattices/Coordinates) and comments.
- **Code Actions**: Automated "Quick Fix" suggestions for INCAR tag typos using Levenshtein distance.
- **Document Symbols**: Comprehensive outline view for POSCAR/CONTCAR and INCAR files.
- **Hover Documentation**: Integrated VASP Wiki documentation for INCAR tags.
- **Auto-completion**: Snippets for KPOINTS and tag suggestions for INCAR.
- **Cross-file Linting**: POTCAR species order validation against POSCAR.
- **GitHub Actions**: Continuous Integration pipeline for linting, building, and testing.

### Fixed

- Robust file URI to path conversion for cross-file validation.
- Resolution of `@typescript-eslint/no-unused-vars` and prettier warnings.
- Duplicate key errors in `package.json`.
- Missing `build` script in `package.json`.

### Changed

- Refactored entry point to use a robust `LspServer` class-based architecture.
- Updated documentation with specialized guides for Neovim setup.

---
*Generated for the VASP Language Server Launch*
