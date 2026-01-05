# DFT Language Server (VASP + CRYSTAL23)

[![Version](https://img.shields.io/badge/version-1.1.0-blue.svg?style=flat-square)](package.json)
[![License](https://img.shields.io/badge/license-MIT-green.svg?style=flat-square)](LICENSE)
[![Support](https://img.shields.io/badge/Target-VASP%206.5.x-EDA523?style=flat-square)](https://www.vasp.at/wiki)

**A strict, high‑performance Language Server for computational materials input files.**

VASP Language Server (a.k.a. the DFT Language Server) brings modern IDE intelligence to:
- **VASP**: `INCAR`, `POSCAR`, `CONTCAR`, `KPOINTS`, `POTCAR`
- **CRYSTAL23**: `.d12`

---

## ⚡ Quick Start

### 1. Install

```bash
npm install -g dft-language-server
```

This installs the `vasp-lsp` and `dft-lsp` CLI entrypoints.

### 2. Connect Your Editor (Neovim)

Add this to your `init.lua` (or use our [Detailed Neovim Guide](docs/NEOVIM.md)):

```lua
require'lspconfig'.vasp_ls.setup{
  cmd = { "vasp-lsp", "--stdio" },
  filetypes = { "vasp", "crystal" },
}
```

**That's it.** Open an `INCAR` file or a `.d12` and enjoy immediate validation.

---

## 🚀 Features

### VASP Support

| Feature | INCAR | POSCAR/CONTCAR | KPOINTS | POTCAR |
| :--- | :---: | :---: | :---: | :---: |
| **Validation** | ✅ | ✅ | ✅ | ✅ |
| **Hover Docs** | ✅ | ✅ | ✅ | - |
| **Completion** | ✅ | - | ✅ | - |
| **Formatting** | ✅ | ✅ | ✅ | - |
| **Quick Fixes** | ✅ | - | - | - |
| **Symbols (Outline)** | ✅ | ✅ | - | - |
| **Semantic Tokens** | ✅ | ✅ | ✅ | - |
| **Folding** | ✅ | ✅ | - | - |

### CRYSTAL23 Support

| Feature | `.d12` |
| :--- | :---: |
| **Validation** | ✅ |
| **Hover Docs** | ✅ |
| **Completion** | ✅ |
| **Semantic Tokens** | ✅ |

### Advanced VASP Support

See [docs/ADVANCED_CONFIGURATION.md](docs/ADVANCED_CONFIGURATION.md) for deep dives into:
- **Scientific precision**: `1E-5`, `.5`, `10*0.0`, Fortran `D` notation
- **Expert workflows**: LDA+U, hybrid functionals, van der Waals
- **Cross‑file intelligence**: POTCAR order checks vs POSCAR/CONTCAR
- **Filetype overrides** for custom naming conventions

---

## 🔧 Usage

### VS Code

Use any generic LSP client extension and point it at the server binary:

- **Command**: `vasp-lsp`
- **Args**: `["--stdio"]`

A minimal VS Code client scaffold is included in `client/vscode` (optional).

### Neovim

We strongly recommend `nvim-lspconfig`. See [docs/NEOVIM.md](docs/NEOVIM.md) for a streamlined setup.

### Other Editors

Any LSP client that can start a stdio server works (Zed, Emacs, Helix, etc.).

---

## 🧰 Custom LSP Commands

The server exposes `workspace/executeCommand` for advanced integrations.

### `vasp.previewStructure`

Return a JSON representation of the parsed structure for the given document.

**Arguments:**
- `{ uri: "file:///path/to/POSCAR" }` (or a plain URI string)

**Result:**
- `{ uri, fileType, structure }` where `structure` is the parsed AST cached by the server

---

## ⚙️ Configuration

The server accepts filetype overrides via LSP settings:

```json
{
  "vasp": {
    "fileTypeOverrides": {
      "filenames": {
        "INCAR.relax": "incar",
        "POSCAR.start": "poscar"
      },
      "extensions": {
        ".vasp": "poscar",
        "kp": "kpoints"
      }
    }
  }
}
```

You can also pass `initializationOptions.fileTypeOverrides` in clients that do not support settings.

---

## 🏗 Architecture & Contribution

We are building a professional-grade tool for the scientific community.

- **Typed, test‑driven**: TypeScript + Jest
- **Performance‑focused**: Robust parsing for 100k+ atom POSCARs
- **Documentation‑driven**: Cross‑validated against the official VASP Wiki

Start here:
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)
- [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md)

---

## 📌 Dependency Notes

Runtime dependencies are purpose-driven:

- **VASP Wiki scraping**: `axios`, `cheerio`, `mwn`
- **HTML → Markdown**: `turndown`, `@joplin/turndown-plugin-gfm`
- **Math rendering**: `mathjax-full`
- **LSP protocol**: `vscode-languageserver`, `vscode-languageserver-textdocument`, `vscode-languageserver-types`

---

## 📄 License

MIT © [TheFermiSea](https://github.com/TheFermiSea)
