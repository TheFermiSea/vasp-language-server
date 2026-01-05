# Neovim Plugin Release Checklist

This repo ships a lightweight Neovim plugin via the `plugin/` and `lua/` directories.
Releases are Git‑based (no package registry).

## Pre‑Release

1. Verify entrypoints exist:
   - `plugin/vasp-lsp.lua`
   - `lua/vasp-lsp.lua`
2. Ensure docs are up to date:
   - `docs/NEOVIM.md`
   - `README.md`
3. Run a quick sanity check in Neovim:
   - `:lua require('vasp-lsp').setup()`
   - Open `INCAR` and `.d12` and confirm server starts
4. Run server tests:
   - `npm test`

## Release

1. Update version/changelog as needed:
   - `package.json` and `CHANGELOG.md`
2. Tag a release:
   - `git tag vX.Y.Z`
3. Push tags:
   - `git push origin vX.Y.Z`

GitHub Actions will:
- Build + test the server
- Create a GitHub release
- Attach the npm tarball (`npm pack` output)

## Post‑Release

- Verify the GitHub Release artifacts
- Announce in release notes
- Update external docs referencing the tag/version
