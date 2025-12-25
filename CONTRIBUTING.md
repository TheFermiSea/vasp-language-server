# Contributing to VASP Language Server

Thank you for your interest in contributing! This project aims to make editing VASP files easier for everyone.

## Development Setup

1. **Prerequisites**: Node.js and npm.
2. **Install**:

    ```bash
    npm install
    ```

3. **Build**:

    ```bash
    npm run build
    ```

    Or watch for changes:

    ```bash
    npm run watch
    ```

## Project Structure

- `src/`: TypeScript source code.
- `src/poscar-parsing.ts`: Where the file structure definition lives.
- `src/poscar-linting.ts`: Where the validation rules live.

## How to Add a New Linting Rule

If you find a common mistake in POSCAR files that isn't caught, here is how to add a check:

1. Open `src/poscar-linting.ts`.
2. Identify the block type you want to check (e.g., `lattice`, `positions`).
3. Find the corresponding function in `poscarBlockLinters`.
4. Add your logic. Use `createDiagnostic` to return an error.

**Example**:

```typescript
// Inside specific linter function
if (someValue < 0) {
    diagnostics.push(createDiagnostic(
        "Value must be positive",
        token.range,
        DiagnosticSeverity.Error
    ));
}
```

## How to Debug

Since the server runs as a subprocess of your editor, debugging can be tricky.

### Method 1: Console Logging

The `logger.ts` wrapper allows you to send logs to your editor's "Output" or "LSP Log" panel.

```typescript
import { logger } from './logger';
logger.info("Checking value: " + value);
```

### Method 2: Node Inspector

You can launch the server with the `--inspect` flag (if your client config allows passing arguments) and attach a generic Node debugger (like Chrome DevTools or VS Code debugger).

## Pull Requests

- Ensure `npm run build` passes without errors.
- Please add comments explaining *why* a change was made if it's complex.
- Update `README.md` if you added a user-facing feature.
