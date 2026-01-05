# DFT Language Server Development Guide

This guide describes the current development workflow for the DFT Language Server (VASP + CRYSTAL23). It supersedes the initial migration notes.

## 1. Project Overview

**Goal**: Provide a fast, strict LSP server for scientific input files:
- **VASP**: `INCAR`, `POSCAR`, `CONTCAR`, `KPOINTS`, `POTCAR`
- **CRYSTAL23**: `.d12`

**Key design goals**:
- Deterministic diagnostics
- High‑performance parsing (100k+ atoms)
- Feature‑based architecture
- Strong test coverage

## 2. Workspace Setup

```bash
npm install
npm run build
```

Run the server locally:

```bash
node out/server.js --stdio
```

Or use the CLI wrapper after building:

```bash
./bin/vasp-lsp --stdio
```

## 3. Scripts & Tooling

- **Build**: `npm run build`
- **Watch**: `npm run watch`
- **Tests**: `npm test`
- **Lint**: `npm run lint`
- **Format**: `npm run format` / `npm run format:check`
- **Update tags**: `npm run update-tags` (network required)

## 4. Repository Layout

```text
src/
  core/            LSP server, cache, parser utils
  features/        File‑type implementations
  data/            Static tag data
  tools/           Standalone scripts (VASP wiki scraper)
  utils/           Shared helpers
  __tests__/       Unit + integration tests
bin/               CLI wrappers
client/            Optional VS Code scaffold
lua/ + plugin/      Neovim plugin
```

## 5. Development Workflow

1. Implement or modify a feature in `src/features/<type>`.
2. Reuse `DocumentCache` where possible to avoid re‑parsing.
3. Wire new providers in `src/core/lsp-server.ts`.
4. Add tests in `src/__tests__/unit` or `src/__tests__/integration`.
5. Run `npm test` and `npm run lint`.

## 6. Debugging

Attach a debugger to the running server using `--inspect`:

```bash
node --inspect=6009 out/server.js --stdio
```

Then attach with a Node debugger (VS Code or other).

## 7. Client Integration

- **Neovim**: See `docs/NEOVIM.md`
- **VS Code**: Optional scaffold in `client/vscode`
- **Other Editors**: Any LSP client that can start a stdio server

## 8. Release Workflow

1. Update `package.json` and `CHANGELOG.md`.
2. Tag a release:
   ```bash
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```
3. GitHub Actions will build, test, and attach an npm tarball to the release.

## 9. Notes on Performance

- POSCAR stress tests live in `src/__tests__/unit/stress-tests.test.ts`.
- Avoid extra allocations in hot parsers (POSCAR/KPOINTS).
- Keep diagnostics precise; prefer shared helpers in `src/utils`.
