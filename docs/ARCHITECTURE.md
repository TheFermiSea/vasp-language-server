# Architecture Guide

This document explains the internal design of the VASP Language Server. It is intended for developers who want to understand how the system works or contribute to the codebase.

## High-Level Design

The project follows a standard **Language Server Protocol (LSP)** client-server architecture.

- **Client**: The editor (e.g., Neovim, VS Code). Handles UI (text buffers, rendering squiggles).
- **Server**: This Node.js process. Handles logic (parsing, validating, documentation).
- **Transport**: Communication happens via `stdio` (Standard Input/Output) using JSON-RPC.

## Directory Structure

```
src/
├── server.ts           # Entry point. Handles LSP connection and events.
├── logger.ts           # Logging wrapper (sends logs to client console).
├── poscar-parsing.ts   # Core Logic: Tokenizes and structures POSCAR files.
├── poscar-linting.ts   # Core Logic: Validates the parsed structural POSCAR data.
├── incar-parsing.ts    # Core Logic: Tokenizes and parses INCAR files (Tags/Values).
├── incar-linting.ts    # Core Logic: Validates INCAR tags against schema.
├── incar-tag.ts        # Data Model: Represents an INCAR tag.
├── potcar-parsing.ts   # Core Logic: Extracts species from POTCAR headers.
├── potcar-linting.ts   # Core Logic: Validates POTCAR against POSCAR.
├── kpoints-parsing.ts  # Core Logic: Parses K-Point generation modes/grids.
├── kpoints-linting.ts  # Core Logic: Validates KPOINTS integrity.
├── structure.ts        # Data Model: Shared interfaces for VASP structures.
├── util.ts             # Generic helpers (string checking, etc.).
└── data/
    └── vasp-tags.ts    # Static database of VASP 6.5.x tag definitions.
```

## Core Components

### 1. The Server (`server.ts`)

- **Role**: The coordinator.
- **Initialization**: Sets up the `createConnection` (from `vscode-languageserver`) and the `TextDocuments` manager. Registers capabilities (Hover, Completion).
- **Event Loop**: Listens for `onDidChangeContent`. When a file changes, it checks the filename (`POSCAR` vs `INCAR`) and delegates to the appropriate validator.
- **Diagnostics**: Collects errors from the Linter and uses `connection.sendDiagnostics` to push them to the client.

### 2. The Parsers

- **POSCAR (`poscar-parsing.ts`)**:
  - **Philosophy**: Fixed line-by-line format.
  - **Mechanism**: State machine expecting specific blocks (Lattice -> Positions).
- **INCAR (`incar-parsing.ts`)**:
  - **Philosophy**: Free-format key-value pairs (`TAG = VALUE`).
  - **Mechanism**: Tokenizer handles semicolons, comments (`#`, `!`), and line continuations (`\`). Groups tokens into logical statements.
- **POTCAR (`potcar-parsing.ts`)**:
  - Extracts `VRHFIN` or `TITEL` headers to identify atomic species.
- **KPOINTS (`kpoints-parsing.ts`)**:
  - Parses generation mode (Line 3) and grid dimensions (Line 4).

### 3. The Linters

- **POSCAR (`poscar-linting.ts`)**:
  - Validates structural constraints (e.g. 3 lattice vectors, correct number of atom coordinates).
- **INCAR (`incar-linting.ts`)**:
  - Validates against `src/data/vasp-tags.ts`.
  - Checks types (Int, Float, Bool, String) and valid string options.
  - Refuses unknown tags by default (warns).
- **POTCAR (`potcar-linting.ts`)**:
  - Cross-references extracted species with `POSCAR` (if present) to detect order mismatches.
- **KPOINTS (`kpoints-linting.ts`)**:
  - Validates that grid dimensions are strictly integers.
  - Checks for recognized generation modes (Monkhorst-Pack, Gamma, etc.).

### 4. Utilities (`util.ts`)

- Pure functions for validating numbers, integers, and strings. Kept simple to be testable and reusable.

## Data Flow (Validation)

1. **User types in Editor**: Client sends `textDocument/didChange`.
2. **Server**: `documents.onDidChangeContent` fires.
3. **Server**: Detects file type.
4. **Parsing**: Delegates to `parsePoscar` or `parseIncar`.
5. **Linting**: Delegates to `validatePoscar` or `validateIncar`.
6. **Response**: Server sends `textDocument/publishDiagnostics` notification.
7. **Client**: Renders warning/error squiggles.
