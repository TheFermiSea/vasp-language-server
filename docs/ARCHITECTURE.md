# Architecture Guide

This document explains the internal design of the VASP Language Server.

## High-Level Design

The project follows a standard **Language Server Protocol (LSP)** client-server architecture.

- **Client**: The editor (e.g., Neovim, VS Code). Handles UI and text rendering.
- **Server**: This Node.js process. Handles parsing, validation, and metadata.
- **Transport**: Standard JSON-RPC over `stdio`.

## Directory Structure

We recently moved to a feature-based modular structure to keep the codebase maintainable:

```bash
src/
├── server.ts           # Entry point. Directs events to feature managers.
├── features/           # LSP Feature Implementations
│   ├── folding.ts      # Folding Range logic
│   ├── incar/          # INCAR-specific features
│   │   ├── completion.ts
│   │   ├── hover.ts
│   │   ├── quick-fix.ts
│   │   └── semantic-tokens.ts
│   └── ... 
├── document-symbols.ts # Outline logic for INCAR/POSCAR.
├── incar-parsing.ts    # INCAR-specific tokenization.
├── poscar-parsing.ts   # POSCAR-specific state machine parser.
├── data/
│   └── vasp-tags.ts    # Static VASP tag database.
└── util.ts             # Shared utilities (Levenstein distance, isNumber, etc.)
```

## Core Components

### 1. The Server (`src/server.ts`)

The server acts as a thin coordinator. It manages the `vscode-languageserver` connection and document synchronization. When an LSP request arrives (e.g., `textDocument/hover`), it identifies the file type and calls the corresponding feature module.

### 2. Feature Modules

Instead of a monolithic linter, features are encapsulated:

- **Semantic Tokens**: Implemented in `src/features/incar/semantic-tokens.ts` using the `SemanticTokensBuilder`.
- **Code Actions**: Typos are caught using Levenshtein distance in `src/features/incar/quick-fix.ts`.
- **Folding**: The `src/features/folding.ts` module parses the file on-the-fly to find logical blocks (like POSCAR coordinate lists).

### 3. Shared Logic

- **Parsers**: Hand-written parsers in `src/*-parsing.ts` convert raw text into structured objects (AST-like) that other features consume.
- **Tag Database**: `src/data/vasp-tags.ts` is the central source of truth for VASP command definitions.

## Data Flow

1. **Change**: Client triggers `didChange`.
2. **Sync**: Server's `TextDocuments` manager updates internal buffer.
3. **Dispatch**: Server calls a validator (e.g., `validateIncar`).
4. **Publish**: Diagnostics are sent back to the client immediately.
5. **Feature Requests**: On-demand requests (Hover, Symbols, etc.) are routed to the relevant module in `src/features/`.

---

### Last Updated: December 2025
