# Neovim Plugin Release Checklist

This repo ships a lightweight Neovim plugin via the `plugin/` and `lua/` directories.
Releases are Git-based (no package registry).

## Pre-Release

1. Verify the entrypoint exists:
   - `plugin/vasp-lsp.lua`
   - `lua/vasp-lsp.lua`
2. Ensure docs are up to date:
   - `docs/NEOVIM.md` includes install + setup instructions.
3. Run a quick sanity check in Neovim:
   - `:lua require('vasp-lsp').setup()`
   - Open `INCAR` and confirm server starts.

## Release

1. Commit changes in the `vasp-language-server` repo.
2. Tag a release:
   - `git tag vX.Y.Z`
3. Push tags:
   - `git push origin vX.Y.Z`

## Post-Release

- Announce in release notes (GitHub Releases).
- Update any external docs referencing the tag/version.
