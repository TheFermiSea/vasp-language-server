# VASP Language Server

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg?style=flat-square)](package.json)
[![License](https://img.shields.io/badge/license-MIT-green.svg?style=flat-square)](LICENSE)
[![Support](https://img.shields.io/badge/Support-VASP%206.5.x-EDA523?style=flat-square)](https://www.vasp.at/wiki)

Intelligent editing experience for Vienna Ab initio Simulation Package (VASP) input files.

[Overview](#overview) • [Features](#features) • [Installation](#installation) • [Usage](#usage) • [Architecture](#architecture) • [Contributing](#contributing)

---

## Overview

The **VASP Language Server** provides modern IDE features to the world of computational materials science. Built on the **Language Server Protocol (LSP)**, it provides instant feedback, validation, and documentation for VASP input files (`INCAR`, `POSCAR`, `KPOINTS`, `POTCAR`).

Stop checking `OUTCAR` for syntax errors. Catch them as you type.

---

## Features

### Feature Matrix

| Feature | INCAR | POSCAR | POTCAR | KPOINTS |
| :--- | :---: | :---: | :---: | :---: |
| **Diagnostics** | ✅ | ✅ | ✅ | ✅ |
| **Hover Docs** | ✅ | - | - | - |
| **Autocomplete** | ✅ | - | - | - |
| **Quick Fixes** | ✅ | - | - | - |
| **Outline (Symbols)** | ✅ | ✅ | - | - |
| **Semantic Tokens** | ✅ | - | - | - |
| **Folding Ranges** | - | ✅ | - | - |

### 🚀 INCAR

* **Smart Validation**: Validates tags against the VASP 6.5.x database.
* **Smart Hover**: View tag descriptions, types, and defaults.
* **Code Actions**: Quick fixes for common typos (e.g., `EMCUT` -> `ENCUT`).
* **Semantic Highlighting**: Rich coloring for tags, numbers, booleans, and comments.

### 📐 POSCAR

* **Lattice Validation**: Checks for 3x3 lattice vectors.
* **Coordinate Check**: Ensures count matches species list.
* **Outline**: Structured view of the file (Lattice, Species, Coordinates).
* **Folding**: Collapse lattice vectors or long coordinate lists.

### 🧪 POTCAR & KPOINTS

* **POTCAR Consistency**: Warns if POTCAR order doesn't match POSCAR.
* **KPOINTS Mode**: Validates generation modes (Monkhorst-Pack, Gamma, etc.).

---

## Installation

### Prerequisites

* **Node.js** (v18 or higher)
* **npm**

### Setup from Source

```bash
git clone https://github.com/TheFermiSea/vasp-language-server.git
cd vasp-language-server
npm install
npm run build
```

---

## Usage

### Running the Server

The server runs via `stdio`:

```bash
node /path/to/repo/out/server.js --stdio
```

### Editor Integration

#### Neovim

We provide a detailed guide for Neovim users (built-in LSP). See [docs/NEOVIM.md](docs/NEOVIM.md).

#### VS Code / Others

Use any generic LSP client extension. Point the command to the path above.

---

## Architecture

The project is modularized by feature. `server.ts` acts as a coordinator delegating to modules in `src/features/`.

* **`src/features/incar/`**: Hover, Completion, Semantic Tokens, Code Actions.
* **`src/features/folding.ts`**: Folding logic.
* **`src/document-symbols.ts`**: Document symbols (Outline).

For more details, see [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

---

## Contributing

We welcome contributions! Please read [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) for development guidelines.

## License

MIT © [TheFermiSea](https://github.com/TheFermiSea)
