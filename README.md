# VASP Language Server

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg?style=flat-square)](package.json)
[![License](https://img.shields.io/badge/license-MIT-green.svg?style=flat-square)](LICENSE)
[![Support](https://img.shields.io/badge/Target-VASP%206.5.x-EDA523?style=flat-square)](https://www.vasp.at/wiki)

**The strictly typed, high-performance Language Server for Computational Materials Science.**

Eliminate `OUTCAR` syntax errors before they happen. VASP Language Server brings modern IDE intelligence to your `INCAR`, `POSCAR`, `KPOINTS`, and `POTCAR` files.

---

## ⚡ Quick Start

### 1. Install Globally

```bash
npm install -g vasp-language-server
```

### 2. Connect Your Editor (Neovim)

Add this to your `init.lua` (or use our [Detailed Neovim Guide](docs/NEOVIM.md)):

```lua
require'lspconfig'.vasp_ls.setup{
  cmd = { "vasp-lsp", "--stdio" },
  filetypes = { "vasp" } -- Ensure you have filetype detection!
}
```

**That's it.** Open an `INCAR` file and enjoy immediate validation.

---

## 🚀 Features

### Core Capabilities

| Feature | INCAR | POSCAR | POTCAR | KPOINTS |
| :--- | :---: | :---: | :---: | :---: |
| **Strict Validation** | ✅ | ✅ | ✅ | ✅ |
| **Hover Documentation** | ✅ | ✅ | - | ✅ |
| **Autocomplete** | ✅ | - | - | ✅ |
| **Quick Fixes** | ✅ | - | - | - |
| **Symbol Outline** | ✅ | ✅ | - | - |
| **Semantic Tokens** | ✅ | ✅ | - | ✅ |
| **Folding** | ✅ | ✅ | - | - |

### Advanced VASP Support <sup>(New in v1.0.0)</sup>

Check out the [Advanced Configuration Guide](docs/ADVANCED_CONFIGURATION.md) for deep dives into:

* **Scientific Precision**: Native support for VASP's flexible numeric formats (e.g., `1E-5`, `.5`, `10*0.0`).
* **Expert Workflows**:
  * **Hubbard U (LDA+U)**: Validates `LDAUL`, `LDAUU`, `LDAUJ` arrays.
  * **Hybrid Functionals**: Full support for `HSE06`, `PBE0` with `HFSCREEN`/`AEXX`.
  * **Van der Waals**: Validation for `IVDW` tags and DFT-D3/D4 methods.
* **Cross-File Intelligence**: Checks `POTCAR` order against `POSCAR` species to prevent severe calculation crashes.

---

## 📦 Installation

**Prerequisites**: Node.js v18 or newer.

```bash
# Install via npm (Recommended)
npm install -g vasp-language-server

# Build from Source
git clone https://github.com/TheFermiSea/vasp-language-server.git
cd vasp-language-server
npm install && npm run build
```

---

## 🔧 Usage

### VS Code

Use any generic LSP client extension (e.g., "Language Server Protocol (LSP) Client") and point it to the binary:

* **Command**: `vasp-lsp`
* **Args**: `["--stdio"]`

### Neovim

We strongly recommend using `nvim-lspconfig`. See [docs/NEOVIM.md](docs/NEOVIM.md) for a zero-friction setup with **LazyVim**.

---

## 🏗 Architecture & Contribution

We are building a professional-grade tool for the scientific community.

* **Written in TypeScript**: Strictly typed for reliability.
* **Zero-Copy Parsing**: Optimized for large `POSCAR` files (100k+ atoms).
* **Documentation-Driven**: Everything is cross-validated against the official [VASP Wiki](https://www.vasp.at/wiki).

See [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) to join us.

---

## 📄 License

MIT © [TheFermiSea](https://github.com/TheFermiSea)
