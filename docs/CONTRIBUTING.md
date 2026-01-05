# Contributing to VASP Language Server

Thank you for your interest in contributing to the VASP Language Server! We aim to provide the best possible developer experience for computational physicists.

## Development Setup

1. **Prerequisites**: Node.js (v18+) and npm.
2. **Install**: `npm install`
3. **Build**: `npm run build`
4. **Test**: `npm test`
5. **Lint & Format**:
   - `npm run lint`: Run ESLint checks (strict no-any policy).
   - `npm run format`: Apply Prettier formatting.
6. **Update Data**:
   - `npm run update-tags`: Fetch latest VASP tags from the official wiki.

## Updating VASP Tags (VASP Wiki)

The `npm run update-tags` command runs the VASP Wiki scraper in `src/tools/vasp-wiki.ts`.
It fetches the latest INCAR tag definitions and writes them to `src/data/vasp-tags.ts`.

**When to run:**
- Before a new release
- When a new VASP version is published
- When users report missing or outdated tags

**What it does:**
1. Queries the VASP Wiki category pages for INCAR tag entries
2. Parses each tag page for type, default, and description fields
3. Regenerates `src/data/vasp-tags.ts`

**Expected output:**
- Console logs showing connection, page counts, and number of tags fetched
- A refreshed `src/data/vasp-tags.ts` file in git status

**Notes:**
- The script is network dependent and will log errors if the wiki is unavailable.
- Please be mindful of rate limits; avoid running it repeatedly in a tight loop.

## Project Structure

We follow a modular, feature-based organization:

- `src/core/lsp-server.ts`: The main coordinator class (`LspServer`).
- `src/core/document-cache.ts`: Centralized AST caching layer.
- `src/features/`: Feature-specific implementations (e.g., `src/features/incar/parsing.ts`).
- `src/utils/`: Shared utilities (`logger.ts`, `util.ts`).
- `src/tools/`: Standalone scripts (`vasp-wiki.ts`).
- `src/__tests__/`: Comprehensive test suites.

## How to Add a New Feature

1. **Implement Logic**: Create your feature in `src/features/` or a new top-level file if it's a core parser.
2. **Utilize Caching**: Ensure your feature uses the `DocumentCache` to retrieve pre-parsed ASTs rather than re-parsing from scratch.
3. **Register in LspServer**: Wire your function into the `LspServer` class handlers (in `src/lsp-server.ts`).
4. **Add Tests**:
   - Add unit tests in `src/__tests__/unit/`.
   - Add integration tests in `src/__tests__/integration/server.test.ts` if appropriate.

## Pull Requests

- **Quality**: Ensure `npm run lint` and `npm test` pass.
- **Performance**: If modifying parsers, run `npm test src/__tests__/unit/stress-tests.test.ts` to ensure no performance regressions.
- **Documentation**: Update [README.md](../README.md) or relevant docs if your feature changes user-facing behavior.

---

### Last Updated: December 2025
