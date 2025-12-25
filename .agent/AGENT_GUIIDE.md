# VASP Language Server Development Agent Guide

This document provides instructions for AI agents working on the VASP Language Server.

## Project Overview

A Language Server Protocol (LSP) implementation for VASP input files (`INCAR`, `POSCAR`).
Written in **TypeScript** using `vscode-languageserver`.
Targeting **VASP 6.5.x** syntax.

## Key Directories

- `src/`: Source code.
  - `server.ts`: Entry point.
  - `*-parsing.ts`: Parsers.
  - `*-linting.ts`: Linters.
  - `data/`: Static data (tags).
- `test/`: Test suite (currently empty, needs implementation).
- `docs/`: User documentation.

## Development Workflow

1. **Beads**: Use `bd` to track issues. Create new issues for new features.
2. **Build**: `npm run build` to compile TypeScript.
3. **Testing**: Run tests (setup pending).
4. **Commits**: Use semantic commit messages.

## Current State

- INCAR/POSCAR validation implemented.
- Testing infrastructure is missing.
- KPOINTS/POTCAR support is missing.
