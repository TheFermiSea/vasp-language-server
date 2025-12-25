# Contributing to VASP Language Server

Thank you for your interest in contributing!

## Development Setup

1. **Prerequisites**: Node.js (v18+) and npm.
2. **Install**: `npm install`
3. **Build**: `npm run build`
4. **Lint & Format**: We use ESLint and Prettier.
   - `npm run lint`: Check for errors.
   - `npm run format`: Format code.

## Project Structure

We follow a feature-based organization:

- `src/features/`: Contains LSP features (Hover, Completion, etc.) grouped by domain.
- `src/*-parsing.ts`: Core parsers for VASP file formats.
- `src/__tests__/`: Unit and integration tests. Please add tests for new features.

## How to Add a New Feature

1. **Define the logic**: Create a new file in `src/features/`.
2. **Register in Server**: Import your function in `src/server.ts` and register it in the appropriate `connection.on...` handler.
3. **Add Tests**: Create a unit test in `src/__tests__/unit/`.

## Pull Requests

- Ensure `npm run lint` and `npm test` pass.
- Maintain a clear and concise coding style.
- Update `README.md` if the change affects users.

---

### Last Updated: December 2025
