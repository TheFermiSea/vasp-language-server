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

## Project Structure

We follow a modular, feature-based organization:

- `src/lsp-server.ts`: The main coordinator class (`LspServer`).
- `src/features/`: Contains LSP feature implementations (Hover, Completion, Semantic Tokens, etc.).
- `src/document-cache.ts`: Centralized AST caching layer.
- `src/*-parsing.ts`: Core parsers for VASP file formats (INCAR, POSCAR, KPOINTS).
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
