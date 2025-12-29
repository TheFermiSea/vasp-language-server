# Advanced Configuration & Features

VASP Language Server (v1.0.0) introduces authoritative validation for complex functional setups. This guide explains the strict scientific rules enforced by the server.

## 1. Hubbard U (LDA+U) Support

The server automatically validates the complex array syntax required for DFT+U calculations.

### Validated Tags

* `LDAU`: Boolean switch.
* `LDAUTYPE`: Integer (1, 2, or 4).
* `LDAUL`: Array of angular moments (l-quantum number).
* `LDAUU`: Array of U values (eV).
* `LDAUJ`: Array of J values (eV).

### Validation Logic

* **Array Length**: Must match `NTYP` (number of species in POSCAR).
* **Type Safety**: `LDAUL` values must be integers (e.g., `-1`, `1`, `2`, `3`).
* **Consistency**: Warns if `LDAU=.TRUE.` but `LDAUU/L` are missing.

## 2. Hybrid Functionals (HSE06 / PBE0)

Hybrid functional setups are notorious for typos. The server strictly checks:

* `LHFCALC = .TRUE.`
* `AEXX`: Exchange fraction (e.g., `0.25` for PBE0).
* `HFSCREEN`: Screening parameter (required for HSE06).
* `AGGAX`, `AGGAC`, `ALDAC`: Advanced exchange-correlation parameters.

## 3. Van der Waals Functionals

Supported `IVDW` (integer) values:

* `1`, `10`: DFT-D2
* `11`, `12`: DFT-D3 (Zero Damping / BJ Damping)
* `4`: dDsC dispersion
* `2`, `20`, `21`: Tkatchenko-Scheffler methods.

## 4. Scientific Number Parsing

VASP is unique in its number parsing. We support:

* **No leading zero**: `.5` is valid.
* **Trailing decimal**: `500.` is valid.
* **Scientific**: `1.0E-5`, `-1.E+02`.
* **Fortran D-notation**: `1.0D-5` (treated as `E`).
* **Repetition**: `10*0.0` (expands to ten zeros).

## 5. Strict Mode (Future)

Currently, all strict validations are **always on**. In future versions, we may introduce a `.vasprc` or client-side setting to relax these checks for legacy files.
