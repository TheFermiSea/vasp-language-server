import {
    createConnection,
    TextDocuments,
    ProposedFeatures,
    InitializeParams,
    TextDocumentSyncKind,
    InitializeResult,
    Diagnostic,
    CompletionItem,
    CompletionItemKind
} from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { validatePoscar } from './poscar-linting';
import { validateIncar } from './incar-linting';
import { validatePotcar } from './potcar-linting';
import { parseKpoints } from './kpoints-parsing';
import { validateKpoints } from './kpoints-linting';
import { parseIncar } from './incar-parsing';
import { formatIncar } from './incar-formatting';
import { logger } from './logger';
import { getIncarCodeActions } from './code-actions';
import { getIncarSymbols, getPoscarSymbols } from './document-symbols';
import { getIncarHover } from './features/incar/hover';
import { getIncarCompletions, resolveIncarCompletion } from './features/incar/completion';
import { getIncarSemanticTokens } from './features/incar/semantic-tokens';
import { semanticTokensLegend } from './features/semantic-tokens-legend';
import { getFoldingRanges } from './features/folding';

/**
 * Create a connection for the server. The connection uses Node's IPC as a transport.
 * Also include all preview / proposed LSP features.
 *
 * @remarks
 * Use `ProposedFeatures.all` to ensure access to the latest LSP capabilities,
 * although we primarily use standard `textDocumentSync` and diagnostics for now.
 */
const connection = createConnection(ProposedFeatures.all);

// Initialize the custom logger wrapper with the actual LSP connection
// so it can send log messages to the client's output channel.
logger.initialize(connection);

/**
 * Create a simple text document manager.
 * The text document manager supports full document sync only.
 *
 * @remarks
 * `TextDocuments` automatically handles `textDocument/didOpen`, `textDocument/onDidChange`,
 * and `textDocument/didClose` notifications, maintaining the state of documents
 * in memory.
 */
const documents = new TextDocuments(TextDocument);

/**
 * Handler for the `initialize` request.
 * This is the first request sent by the client (VS Code, Neovim, etc.) to the server.
 *
 * @param params - Initialization parameters containing client capabilities.
 * @returns An `InitializeResult` containing the server's capabilities.
 */
connection.onInitialize((_params: InitializeParams) => {
    logger.info('Initializing VASP Language Server...');
    logger.info('[Info] - Mic check... one, two... Tuning lattice vectors to A=440Hz...');
    logger.info('Happy Git-Mitzvah @TraceBivens! You are finally a post-pubescent coder!');

    // Define the capabilities the server supports
    const result: InitializeResult = {
        capabilities: {
            // Incremental sync is more efficient than full sync for large files
            textDocumentSync: TextDocumentSyncKind.Incremental,

            // Hover support (e.g. documentation for tags)
            hoverProvider: true,

            // Completion support
            completionProvider: {
                resolveProvider: true,
                triggerCharacters: ['=']
            },

            // Formatting support
            documentFormattingProvider: true,

            // Code Action support (Quick Fixes)
            codeActionProvider: true,

            // Document Symbol support (Outline)
            documentSymbolProvider: true,

            // Semantic Tokens (Syntax Highlighting)
            semanticTokensProvider: {
                legend: semanticTokensLegend,
                full: true
            },

            // Folding Ranges
            foldingRangeProvider: true
        }
    };
    return result;
});

/**
 * Handler for the `initialized` notification.
 * This is sent after the client has received the result of the `initialize` request.
 */
connection.onInitialized(() => {
    logger.info('VASP Language Server Initialized.');
});

/**
 * Validation Listener
 * Triggered whenever a text document is changed (or opened).
 * This is the entry point for the linting workflow.
 */
documents.onDidChangeContent((change) => {
    validateTextDocument(change.document);
});

/**
 * Validates a text document by parsing it and checking for errors.
 *
 * @param textDocument - The document to validate.
 * @returns A promise that resolves when validation is complete.
 */
async function validateTextDocument(textDocument: TextDocument): Promise<void> {
    const uri = textDocument.uri;
    const diagnostics: Diagnostic[] = [];

    const fileName = uri.split('/').pop() || '';

    // Simple heuristic to check file type
    if (fileName.match(/POSCAR/i) || fileName.match(/CONTCAR/i)) {
        try {
            logger.info(`Validating POSCAR: ${uri}`);
            const poscarDiagnostics = validatePoscar(textDocument);
            diagnostics.push(...poscarDiagnostics);
        } catch (e) {
            logger.error(`Error validating POSCAR: ${e}`);
        }
    } else if (fileName.match(/INCAR/i)) {
        try {
            logger.info(`Validating INCAR: ${uri}`);
            // Parse and Lint INCAR
            const parsed = parseIncar(textDocument);
            const incarDiagnostics = validateIncar(parsed);
            diagnostics.push(...incarDiagnostics);
        } catch (e) {
            logger.error(`Error validating INCAR: ${e}`);
        }
    } else if (fileName.match(/POTCAR/i)) {
        try {
            logger.info(`Validating POTCAR: ${uri}`);
            const potcarDiagnostics = validatePotcar(textDocument);
            diagnostics.push(...potcarDiagnostics);
        } catch (e) {
            logger.error(`Error validating POTCAR: ${e}`);
        }
    } else if (fileName.match(/KPOINTS/i)) {
        try {
            logger.info(`Validating KPOINTS: ${uri}`);
            const parsed = parseKpoints(textDocument);
            const kpointsDiagnostics = validateKpoints(parsed);
            diagnostics.push(...kpointsDiagnostics);
        } catch (e) {
            logger.error(`Error validating KPOINTS: ${e}`);
        }
    }

    // Send the computed diagnostics to the client.
    // If 'diagnostics' is empty, this clears any previous errors for the file.
    connection.sendDiagnostics({ uri: textDocument.uri, diagnostics });
}

// Support for Completion
connection.onCompletion((_textDocumentPosition: any): CompletionItem[] => {
    const uri = _textDocumentPosition.textDocument.uri;
    const items: CompletionItem[] = [];

    // KPOINTS Snippets
    if (uri.match(/KPOINTS/i)) {
        items.push({
            label: 'Monkhorst-Pack',
            kind: CompletionItemKind.Snippet,
            detail: 'Generate automatic Monkhorst-Pack grid',
            insertText: 'Automatic Mesh\n0\nMonkhorst-Pack\n4 4 4\n0 0 0'
        });
        items.push({
            label: 'Gamma-Centered',
            kind: CompletionItemKind.Snippet,
            detail: 'Generate automatic Gamma-centered grid',
            insertText: 'Automatic Mesh\n0\nGamma\n4 4 4\n0 0 0'
        });
        return items;
    }

    // INCAR Completions (Default)
    return getIncarCompletions();
});

connection.onCompletionResolve((item: CompletionItem): CompletionItem => {
    return resolveIncarCompletion(item);
});

// Support for Hover
connection.onHover(({ textDocument, position }) => {
    const document = documents.get(textDocument.uri);
    if (!document) return null;

    // Check if INCAR
    if (!textDocument.uri.match(/INCAR/i)) return null;

    return getIncarHover(document, position);
});

// Support for Formatting
connection.onDocumentFormatting((params) => {
    const document = documents.get(params.textDocument.uri);
    if (!document) return [];

    // Only format INCAR for now
    if (document.uri.match(/INCAR/i)) {
        return formatIncar(document);
    }
    return [];
});

// Support for Code Actions
connection.onCodeAction((params) => {
    return getIncarCodeActions(params.textDocument.uri, params.context.diagnostics);
});

// Support for Semantic Tokens
connection.languages.semanticTokens.on((params) => {
    const document = documents.get(params.textDocument.uri);
    if (!document) return { data: [] };

    if (document.uri.match(/INCAR/i)) {
        return getIncarSemanticTokens(document);
    }
    return { data: [] };
});

// Support for Document Symbols (Outline)
connection.onDocumentSymbol((params) => {
    const document = documents.get(params.textDocument.uri);
    if (!document) return null;

    if (document.uri.match(/INCAR/i)) {
        return getIncarSymbols(document);
    } else if (document.uri.match(/POSCAR/i) || document.uri.match(/CONTCAR/i)) {
        return getPoscarSymbols(document);
    }

    return null;
});

// Make the text document manager listen on the connection
documents.listen(connection);

// Support for Folding Ranges
connection.onFoldingRanges((params) => {
    const document = documents.get(params.textDocument.uri);
    if (!document) return [];
    return getFoldingRanges(document);
});

// Listen on the connection
connection.listen();
