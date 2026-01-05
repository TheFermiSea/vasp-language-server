# VASP Language Server (VS Code)

This extension launches the VASP Language Server using the `vasp-lsp` binary.

## Prerequisites

- Install the server globally: `npm install -g vasp-language-server`
- Ensure `vasp-lsp` is on your PATH

## Configuration

You can override the server path in settings:

```json
{
  "vaspLsp.serverPath": "vasp-lsp"
}
```

## File Support

The extension activates on VASP filenames:

- `INCAR`
- `POSCAR`
- `CONTCAR`
- `POTCAR`
- `KPOINTS`
