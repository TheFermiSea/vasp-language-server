# Advanced Configuration & Features

This guide covers advanced VASP validations and server configuration for the DFT Language Server (v1.1.0).

## Supported File Types

- **VASP**: `INCAR`, `POSCAR`, `CONTCAR`, `KPOINTS`, `POTCAR`
- **CRYSTAL23**: `.d12`

## 1. Hubbard U (LDA+U) Support

The server validates the complex array syntax required for DFT+U calculations.

**Validated tags**
- `LDAU`: Boolean switch
- `LDAUTYPE`: Integer (1, 2, or 4)
- `LDAUL`: Array of angular moments (l‑quantum number)
- `LDAUU`: Array of U values (eV)
- `LDAUJ`: Array of J values (eV)

**Validation logic**
- **Array length** must match `NTYP` (number of species in POSCAR).
- **Type safety**: `LDAUL` values must be integers (e.g., `-1`, `1`, `2`, `3`).
- **Consistency**: Warns if `LDAU=.TRUE.` but `LDAUU/L` are missing.

## 2. Hybrid Functionals (HSE06 / PBE0)

Hybrid functional setups are prone to typos. The server checks:

- `LHFCALC = .TRUE.`
- `AEXX`: Exchange fraction (e.g., `0.25` for PBE0)
- `HFSCREEN`: Screening parameter (required for HSE06)
- `AGGAX`, `AGGAC`, `ALDAC`: Advanced exchange‑correlation parameters

## 3. Van der Waals Functionals

Supported `IVDW` values:

- `1`, `10`: DFT‑D2
- `11`, `12`: DFT‑D3 (Zero Damping / BJ Damping)
- `4`: dDsC dispersion
- `2`, `20`, `21`: Tkatchenko‑Scheffler methods

## 4. Scientific Number Parsing

VASP accepts several numeric formats. We support:

- **No leading zero**: `.5`
- **Trailing decimal**: `500.`
- **Scientific**: `1.0E-5`, `-1.E+02`
- **Fortran D‑notation**: `1.0D-5` (treated as `E`)
- **Repetition**: `10*0.0` (expands to ten zeros)

## 5. Cross‑File POTCAR Validation

When editing `POTCAR`, the server checks element ordering against `POSCAR`. If `POSCAR` is missing, it will fall back to `CONTCAR`.

## 6. Filetype Overrides

If your files use custom names or extensions, override detection via LSP settings.

**Supported types:** `incar`, `poscar`, `kpoints`, `potcar`, `crystal`

Example (Neovim `lspconfig`):

```lua
lspconfig.vasp_ls.setup{
  settings = {
    vasp = {
      fileTypeOverrides = {
        filenames = {
          ["INCAR.relax"] = "incar",
          ["POSCAR.start"] = "poscar",
        },
        extensions = {
          [".vasp"] = "poscar",
          ["kp"] = "kpoints",
        },
      },
    },
  },
}
```

Notes:
- Filename matches are case‑insensitive.
- Extensions may include or omit the leading dot.

You can also pass overrides via `initializationOptions.fileTypeOverrides` for clients that do not support settings.

## 7. Custom LSP Commands

The server exposes `workspace/executeCommand` for advanced tooling.

### `vasp.previewStructure`

Returns the parsed AST for the given document (useful for debugging and integrations).

**Arguments:**
- `{ uri: "file:///path/to/POSCAR" }` or a raw URI string

**Result:**
- `{ uri, fileType, structure }`

## 8. Strict Mode (Future)

Currently, strict validations are always on. A future release may add a server setting to relax checks for legacy files.
