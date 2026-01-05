# Contributing to the DFT Language Server

Thanks for your interest in contributing! This project targets both **VASP** and **CRYSTAL23** input files with a strict, fast LSP implementation.

## Development Setup

1. **Prerequisites**: Node.js v18+ and npm
2. **Install**: `npm install`
3. **Build**: `npm run build`
4. **Test**: `npm test`
5. **Lint & Format**:
   - `npm run lint`
   - `npm run format` (or `npm run format:check` in CI)

## Updating VASP Tags (VASP Wiki)

`npm run update-tags` runs the VASP Wiki scraper in `src/tools/vasp-wiki.ts`.
It fetches INCAR tags and regenerates `src/data/vasp-tags.ts`.

**When to run:**
- Before a release
- When a new VASP version ships
- When users report missing/outdated tags

**Notes:**
- The script is network‑dependent.
- Be mindful of rate limits when iterating.

## Project Structure

- `src/core/`: LSP server, cache, and utilities
- `src/features/`: File‑type implementations (incar/poscar/kpoints/potcar/crystal)
- `src/data/`: Static tag data
- `src/tools/`: Standalone scripts (e.g., `vasp-wiki.ts`)
- `src/__tests__/`: Unit + integration tests

## Adding a New Feature

1. **Implement logic** in a new or existing `src/features/<type>` module.
2. **Use the cache** where possible to avoid duplicate parsing.
3. **Wire it up** in `src/core/lsp-server.ts` (provider map).
4. **Add tests**:
   - Unit tests in `src/__tests__/unit/`
   - Integration tests in `src/__tests__/integration/` if appropriate

## Pull Requests

- Ensure `npm run lint` and `npm test` pass.
- For parser changes, run the stress tests:
  - `npm test -- src/__tests__/unit/stress-tests.test.ts`
- Update docs if user‑facing behavior changes.

## Release Workflow

Tagging a version triggers the GitHub Actions release workflow:

- `git tag vX.Y.Z`
- `git push origin vX.Y.Z`

The workflow builds, tests, and attaches an npm tarball to the GitHub Release.

---

### Last Updated: January 2026
