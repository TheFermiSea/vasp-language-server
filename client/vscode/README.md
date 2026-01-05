# DFT Language Server (VS Code)

This optional scaffold launches the DFT Language Server using the `vasp-lsp` (or `dft-lsp`) binary.

## Prerequisites

- Install the server globally: `npm install -g dft-language-server`
- Ensure `vasp-lsp` is on your PATH

## Configuration

You can override the server path in settings:

```json
{
  "vaspLsp.serverPath": "vasp-lsp"
}
```

## File Support

The extension activates on common VASP filenames and CRYSTAL23 `.d12`:

- `INCAR`
- `POSCAR`
- `CONTCAR`
- `POTCAR`
- `KPOINTS`
- `*.d12`

## Notes

This is intentionally minimal. Most users should use a generic LSP client or Neovim setup from `docs/NEOVIM.md`.
