# Architecture Guide

This document explains the internal design of the DFT Language Server (VASP + CRYSTAL23).

## High-Level Design

The project follows a standard **Language Server Protocol (LSP)** client‑server architecture, optimized for performance and low memory overhead when handling large scientific datasets (e.g., massive POSCAR files).

- **Client**: The editor (Neovim, VS Code, Zed, Emacs, etc.)
- **Server**: Node.js process (`vasp-lsp` / `dft-lsp`)
- **Transport**: JSON‑RPC over `stdio`

## Directory Structure

We use a feature‑based structure with a centralized caching layer:

```bash
src/
├── core/
│   ├── lsp-server.ts        # LspServer class (main coordinator)
│   ├── feature-provider.ts  # FeatureProvider registry + base class
│   ├── parser-utils.ts      # Shared parsing helpers
│   └── document-cache.ts    # LRU cache for parsed structures
├── features/               # Feature-specific logic
│   ├── incar/              # INCAR features (completion, hover, formatting)
│   ├── poscar/             # POSCAR/CONTCAR features (linting, symbols)
│   ├── kpoints/            # KPOINTS features (completion, formatting)
│   ├── potcar/             # POTCAR features (async linting)
│   └── crystal/            # CRYSTAL23 features (.d12)
├── data/                   # Static tag data (vasp-tags.ts, crystal-tags.ts)
├── utils/                  # Shared utilities (logger, diagnostics, ranges)
├── tools/                  # Standalone tools (vasp-wiki.ts)
└── server.ts               # Entry point
```

Additional packaging/artifacts:

```bash
bin/          # CLI wrappers (vasp-lsp, dft-lsp)
client/       # Optional VS Code scaffold
lua/ + plugin/ # Neovim plugin
```

## Core Components

### 1. The LspServer (`src/core/lsp-server.ts`)

`LspServer` is the central coordinator:

- Initializes the LSP connection and registers handlers
- Detects file types from URIs (including `.d12`)
- Applies filetype overrides (settings + initialization options)
- Routes requests to the correct `FeatureProvider`
- Exposes custom commands such as `vasp.previewStructure`

### 2. FeatureProvider Registry (`src/core/feature-provider.ts`)

Each file type has a provider that implements hover, completion, linting, formatting, etc. The registry avoids switch statements and keeps feature logic modular.

### 3. Document Cache (`src/core/document-cache.ts`)

Parsed structures are cached in an **LRU** cache to avoid repeat work across features (e.g., folding, symbols, and hover). Entries are invalidated on document change.

### 4. Parser Utilities (`src/core/parser-utils.ts`)

Shared helpers for tokenization, ranges, and safe parsing reduce duplication and keep feature logic consistent.

### 5. Asynchronous Linting

Cross‑file validation (e.g., POTCAR vs POSCAR/CONTCAR) is async to avoid blocking the LSP thread.

## Architecture Diagram (Request Flow)

```mermaid
flowchart TD
    A[LSP Client] -->|Request| B[LspServer]
    B --> C{Detect FileType}
    C -->|incar| D[IncarFeatureProvider]
    C -->|poscar| E[PoscarFeatureProvider]
    C -->|kpoints| F[KpointsFeatureProvider]
    C -->|potcar| G[PotcarFeatureProvider]
    C -->|crystal| H[CrystalFeatureProvider]
    D --> I[Parser/Linter/Features]
    E --> I
    F --> I
    G --> I
    H --> I
    B <--> J[DocumentCache (LRU)]
    I <--> J
    I -->|Diagnostics/Edits| K[LSP Response]
```

## Data Flow

1. **didChange**: Client sends an update.
2. **Detect file type**: `LspServer` selects the matching provider (or override).
3. **Cache update**: Providers parse and update `DocumentCache`.
4. **Lint**: Diagnostics are generated (async where needed).
5. **Publish**: Diagnostics are sent to the client.
6. **Feature requests**: Hover, symbols, tokens, etc. read from cache.

---

### Last Updated: January 2026
