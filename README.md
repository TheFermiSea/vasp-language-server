# VASP Language Server

A Language Server Protocol (LSP) implementation for VASP (Vienna Ab initio Simulation Package) input files (`POSCAR`, `INCAR`, `POTCAR`, `KPOINTS`).

![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

## Overview

This project provides intelligent editing features for VASP files in any LSP-compliant editor (VS Code, Neovim, Zed, Emacs, etc.). It supports robust validation for `POSCAR` structural files and advanced `INCAR` editing with type checking and documentation.

## Features

- **INCAR Support (VASP 6.5.x)**:
  - **Validation**: Strict type checking for 50+ common tags (e.g., `ISMEAR` must be Int, `LREAL` must be Boolean).
  - **Hover Documentation**: Detailed descriptions and default values for tags (e.g., valid flags for `ALGO`).
  - **Autocomplete**: Suggestions for known tags and boolean values (.TRUE./.FALSE.).
  - **Syntax Checks**: Detects missing equals signs, invalid quoting, and continuation line issues.
  - **Parsing**: Handles complex syntax including inline comments (`#`, `!`), semicolons (`;`), and line continuations (`\`).

- **POSCAR Validation**:
  - Validates precise file structure (Scaling factor, Lattice vectors, Atom counts).
  - Checks for correct number of tokens per line.
  - Validates constraints (e.g., lattice vectors must be 3 numbers, atom counts must be integers).
  - Cross-validation: Ensures the number of species names matches the number of atom counts.
  - Supports VASP 5.x format (implicit species names) and VASP 4.x.
  - Validates Selective Dynamics ('T'/'F' flags).

## Installation

### Prerequisites

- Node.js (>= 14.0)
- npm

### Building from Source

1. Clone the repository:

    ```bash
    git clone https://github.com/TheFermiSea/vasp-language-server.git
    cd vasp-language-server
    ```

2. Install dependencies:

    ```bash
    npm install
    ```

3. Build the server:

    ```bash
    npm run build
    ```

    This compiles the TypeScript source to `out/server.js`.

4. (Optional) Create the binary wrapper:

    ```bash
    npm pkg set bin.vasp-lsp="./bin/vasp-lsp"
    chmod +x bin/vasp-lsp
    ```

## Usage

### Neovim (nvim-lspconfig)

Add the following configuration to your `init.lua` or LSP config file:

```lua
local lspconfig = require('lspconfig')
local configs = require('lspconfig.configs')

if not configs.vasp_lsp then
  configs.vasp_lsp = {
    default_config = {
      cmd = { 'node', '/path/to/vasp-language-server/out/server.js', '--stdio' },
      filetypes = { 'vasp', 'poscar', 'incar' },
      root_dir = lspconfig.util.root_pattern('.git', 'INCAR', 'POSCAR'),
      settings = {},
    },
  }
end

lspconfig.vasp_lsp.setup{}
```

### VS Code

To use this server with VS Code, you typically need a client extension wrapper.
(See `client/` folder if available, or use a generic LSP client extension like "Run on Save").

## Contributing

See [CONTRIBUTING.md](./CONTRIBUTING.md) for details on how to set up the development environment and add new features.
