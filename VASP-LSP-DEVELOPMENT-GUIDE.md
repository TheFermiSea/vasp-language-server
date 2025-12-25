# VASP LSP Development Guide

## 1. Project Overview & Objective

**Goal**: Create a generic Language Server Protocol (LSP) server for VASP files (specifically `POSCAR` and `INCAR`) by porting logic from the `vscode-vasp-support` VS Code extension.

**Target**: The resulting server must run as a standalone Node.js process, communicable via stdio JSON-RPC, usable by editors like Neovim, Zed, VS Code, or Emacs.

### Features

The LSP should provide the following capabilities:

- **Validation/Diagnostics**: Real-time linting of `POSCAR` and `INCAR` files.
- **Hover Documentation**: Show explanations for INCAR tags (e.g., `ENCUT`, `ISMEAR`) when hovering.
- **Completion** (Future): Auto-completion for INCAR tags.
- **Formatting** (Future): Auto-formatting of input files.

---

## 2. Workspace Setup

Initialize a new directory `vasp-language-server`.

### 2.1. Initialization

Run the following commands to scaffold the project:

```bash
mkdir vasp-language-server
cd vasp-language-server
npm init -y
```

**Install LSP Dependencies:**

```bash
npm install vscode-languageserver vscode-languageserver-textdocument vscode-languageserver-types
```

**Install Development Dependencies:**

```bash
npm install -D typescript @types/node eslint prettier @typescript-eslint/parser @typescript-eslint/eslint-plugin
```

**Install Dependencies required by original logic:**

```bash
npm install axios cheerio turndown
npm install -D @types/turndown
```

### 2.2. Configuration Files

#### `tsconfig.json`

Create a strict configuration to ensure type safety during the port.

```json
{
  "compilerOptions": {
    "module": "commonjs",
    "target": "es2020",
    "outDir": "out",
    "rootDir": "src",
    "lib": ["es2020"],
    "sourceMap": true,
    "strict": true,
    "moduleResolution": "node",
    "esModuleInterop": true,
    "skipLibCheck": true,
    "forceConsistentCasingInFileNames": true
  },
  "include": ["src"]
}
```

#### `.gitignore`

Standard node ignore file.

```text
node_modules/
out/
.DS_Store
```

---

## 3. Migration Strategy (File by File)

You will migrate files from the `vscode-vasp-support` extension to `src/` in the new project.

### Phase A: Direct Copies (With Minor Tweaks)

Copy these files directly. They contain pure logic or standard Node libraries and should require minimal changes.

1. **`src/util.ts`**
    - **Action**: Copy file.
    - **Check**: Ensure no `vscode` imports exist.
2. **`src/incar-tag.ts`**
    - **Action**: Copy file.
    - **Check**: This should be pure interface definitions (safe).
3. **`src/html-to-markdown.ts`**
    - **Action**: Copy file.
4. **`src/vasp-wiki.ts`**
    - **Action**: Copy file.
    - **Modification**: This file likely uses `logger.ts`. Ensure references to `log()` or `error()` match the new logger implementation (Phase B).

### Phase B: Infrastructure Replacements

#### 1. `src/logger.ts`

The original logger likely writes to a `vscode.OutputChannel`. The LSP equivalent is sending log messages to the client via the connection.

**New Implementation:**

```typescript
import { Connection } from 'vscode-languageserver';

let connection: Connection | null = null;

export function setConnection(conn: Connection) {
    connection = conn;
}

export function log(message: string) {
    connection?.console.log(message);
}

export function error(message: string) {
    connection?.console.error(message);
}
```

### Phase C: Core Logic Refactoring (Critical)

These files rely heavily on `vscode` types (like `vscode.Range`, `vscode.Diagnostic`). You must replace them with `vscode-languageserver-types`.

#### 1. `src/poscar-parsing.ts`

- **Imports**: Remove `import * as vscode from 'vscode'`. Add `import { Range, Position } from 'vscode-languageserver-types'`.
- **TextLine Replacement**: `vscode.TextLine` does not exist in LSP server libraries.
  - **Strategy**: Create a simplified interface or pass raw text lines.
  - **Refactor**: Redefine `PoscarLine`.

```typescript
import { Range } from 'vscode-languageserver-types';

export interface PoscarLine {
    text: string;
    lineNumber: number; // 0-based
    range: Range;
}
```

- **Helper**: Use a helper to split full document text into these lines.

```typescript
import { TextDocument } from 'vscode-languageserver-textdocument';

export function getDocumentLines(document: TextDocument): PoscarLine[] {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
    return lines.map((text, i) => ({
        text,
        lineNumber: i,
        range: {
            start: { line: i, character: 0 },
            end: { line: i, character: text.length }
        }
    }));
}
```

#### 2. `src/poscar-linting.ts`

- **Imports**: Remove `vscode`. Import `{ Diagnostic, DiagnosticSeverity, Range }` from `vscode-languageserver-types`.
- **Diagnostic Creation**:
  - Change return type from `vscode.Diagnostic` to `Diagnostic`.
  - **Severity Mapping**:
    - `vscode.DiagnosticSeverity.Error` -> `DiagnosticSeverity.Error`
    - `vscode.DiagnosticSeverity.Warning` -> `DiagnosticSeverity.Warning`
- **Validation Entry Point**: Refactor the main function to accept a `TextDocument` (from `vscode-languageserver-textdocument`) instead of `vscode.TextDocument`. Return `Diagnostic[]`.

### Phase D: The Server Entry Point (`src/server.ts`)

Create this **NEW** file to replace `extension.ts`. This wires everything together.

```typescript
import {
    createConnection,
    TextDocuments,
    ProposedFeatures,
    InitializeParams,
    DidChangeConfigurationNotification,
    TextDocumentSyncKind,
    InitializeResult,
    Hover,
    Diagnostic
} from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { validatePoscar } from './poscar-linting'; // Your refactored linter
// import { getIncarHover } from './vasp-wiki'; // If implementing Hover

// Create connection via stdio
const connection = createConnection(ProposedFeatures.all);

// Initialize Logger
import { setConnection } from './logger';
setConnection(connection);

// Text Document Manager
const documents = new TextDocuments(TextDocument);

connection.onInitialize((params: InitializeParams) => {
    const result: InitializeResult = {
        capabilities: {
            textDocumentSync: TextDocumentSyncKind.Incremental,
            hoverProvider: true, // Set to true if implementing hover
        }
    };
    return result;
});

// Validate on change
documents.onDidChangeContent(change => {
    validateTextDocument(change.document);
});

async function validateTextDocument(textDocument: TextDocument): Promise<void> {
    const uri = textDocument.uri;
    const diagnostics: Diagnostic[] = [];

    // Simple file extension check (or rely on client sending only correct files)
    if (uri.match(/POSCAR/i)) {
         // Assume validatePoscar parses and returns Diagnostic[]
         // You might need to catch errors to prevent server crash
         try {
             // const poscarDiagnostics = validatePoscar(textDocument);
             // diagnostics.push(...poscarDiagnostics);
         } catch (e) {
             connection.console.error(`Error validating POSCAR: ${e}`);
         }
    }

    connection.sendDiagnostics({ uri: textDocument.uri, diagnostics });
}

documents.listen(connection);
connection.listen();
```

---

## 4. Testing & Verification

### Building

```bash
npm run build
# OR
npx tsc
```

### Manual Testing with Stdio

You can technically run the server manually and pipe JSON-RPC to it, but it's hard. Better to use a client.

```bash
node out/server.js --stdio
```

### Debugging in VS Code

To debug *this* server code while it runs inside a client (like VS Code extension or even just attached), add a `.vscode/launch.json`:

```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "type": "node",
            "request": "attach",
            "name": "Attach to Server",
            "port": 6009,
            "restart": true,
            "sourceMaps": true
        }
    ]
}
```

*Note: You will need to instruct the *Client* to start the server with `--inspect=6009`.*

---

## 5. Client Integration Guide

### 5.1 Neovim (Lua)

Once built, configure Neovim's built-in LSP client.

`~/.config/nvim/init.lua`:

```lua
local lspconfig = require('lspconfig')
local configs = require('lspconfig.configs')

if not configs.vasp_lsp then
  configs.vasp_lsp = {
    default_config = {
      cmd = { 'node', '/path/to/vasp-language-server/out/server.js', '--stdio' },
      filetypes = { 'vasp', 'poscar', 'incar' },
      root_dir = lspconfig.util.root_pattern('.git', 'INCAR', 'POSCAR'),
      settings = {},
    },
  }
end

lspconfig.vasp_lsp.setup{}
```

### 5.2 VS Code Extension (Client)

To use this in VS Code, you typically create a *separate* minimal extension that starts this server.

1. Create a folder `vasp-client`.
2. Install `vscode-languageclient`.
3. In `extension.ts`:

```typescript
import { LanguageClient, ServerOptions, TransportKind } from 'vscode-languageclient/node';

// ... inside activate() ...
const serverModule = context.asAbsolutePath(path.join('..', 'vasp-language-server', 'out', 'server.js'));

const serverOptions: ServerOptions = {
    run: { module: serverModule, transport: TransportKind.stdio },
    debug: {
        module: serverModule,
        transport: TransportKind.stdio,
        options: { execArgv: ['--nolazy', '--inspect=6009'] }
    }
};

// ... create client and start ...
```

---

## 6. Packaging

To distribute:

1. **Publish to NPM**: `npm publish` (Ensure `bin` entry in package.json points to a CLI wrapper).
2. **Bin Wrapper** (`bin/vasp-lsp`):

    ```javascript
    #!/usr/bin/env node
    require('../out/server.js');
    ```

3. Make executable: `chmod +x bin/vasp-lsp`.
