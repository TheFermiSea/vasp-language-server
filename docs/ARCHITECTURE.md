# Architecture Guide

This document explains the internal design of the VASP Language Server.

## High-Level Design

The project follows a standard **Language Server Protocol (LSP)** client-server architecture, optimized for high performance and low memory overhead when handling large scientific datasets (e.g., massive POSCAR files).

- **Client**: The editor (e.g., Neovim, VS Code).
- **Server**: This Node.js process. Parses VASP files and provides intelligent feedback.
- **Transport**: JSON-RPC over `stdio`.

## Directory Structure

We use a feature-based modular structure with a centralized caching layer:

```bash
src/
├── core/
│   ├── lsp-server.ts       # Core LspServer class (Main coordinator)
│   └── document-cache.ts   # AST Caching layer (Standardized VaspStructure)
├── features/               # Feature-specific Logic
│   ├── incar/              # INCAR features (Completion, Hover, Folding, etc.)
│   ├── poscar/             # POSCAR features (Folding, Hover, Symbols)
│   ├── kpoints/            # KPOINTS features (Snippets, Hover)
│   └── potcar/             # POTCAR features (Async Linting)
├── utils/                  # Shared Utilities
├── tools/                  # Standalone Tools (vasp-wiki.ts)
└── server.ts               # Entry point (Thin wrapper)
```

## Core Components

### 1. The LspServer (`src/lsp-server.ts`)

The `LspServer` class is the central nervous system. It:

- Manages the `_Connection`.
- Orchestrates `validateTextDocument` across different VASP formats.
- Dispatches events to optimized feature managers.
- **Strictly Typed**: Uses `vscode-languageserver` interfaces for all handlers.

### 2. Document Cache (`src/document-cache.ts`)

To prevent redundant parsing across multiple LSP features (e.g., Outline and Folding both needing the same POSCAR AST), we use a `DocumentCache`.

- It stores a `VaspStructure` (a union of `IncarDocument`, `PoscarDocument`, or `KpointsDocument`).
- It is automatically invalidated when file content changes.

### 3. High-Performance Parsers

Our parsers are hand-written for speed and memory efficiency:

- **POSCAR**: Uses a specialized state machine to handle files with 100,000+ atoms in <400ms.
- **INCAR**: Converts assignments into a key-value AST with support for comments and types.

### 4. Asynchronous Linting

Features like **POTCAR/POSCAR consistency checks** are implemented asynchronously (`src/potcar-linting.ts`) to ensure the main LSP thread never blocks on disk I/O, keeping the editor responsive.

## Data Flow

1. **didChange**: Client sends update.
2. **Standardize**: `LspServer` parses the document and updates the `DocumentCache`.
3. **Linter**: `LspServer` runs the relevant linter (e.g., `validateIncar`) using the cached AST.
4. **Publish**: Diagnostics are pushed to the client.
5. **Feature Requests**: On-demand requests (Hover, Symbols, etc.) pull the pre-parsed AST from the cache and return results instantly.

---

### Last Updated: December 2025
