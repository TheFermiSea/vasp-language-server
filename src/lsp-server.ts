import {
    createConnection,
    TextDocuments,
    ProposedFeatures,
    InitializeParams,
    TextDocumentSyncKind,
    InitializeResult,
    Diagnostic,
    CompletionItem,
    CompletionItemKind,
    _Connection
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
 * Encapsulates the VASP Language Server logic.
 */
export class LspServer {
    private connection: _Connection;
    private documents: TextDocuments<TextDocument>;

    constructor() {
        this.connection = createConnection(ProposedFeatures.all);
        this.documents = new TextDocuments(TextDocument);

        // Initialize logger with connection
        logger.initialize(this.connection);

        this.registerHandlers();
    }

    /**
     * Registers all LSP event handlers.
     */
    private registerHandlers(): void {
        this.connection.onInitialize(this.onInitialize.bind(this));
        this.connection.onInitialized(this.onInitialized.bind(this));

        this.documents.onDidChangeContent((change) => {
            this.validateTextDocument(change.document);
        });

        this.connection.onCompletion(this.onCompletion.bind(this));
        this.connection.onCompletionResolve(this.onCompletionResolve.bind(this));
        this.connection.onHover(this.onHover.bind(this));
        this.connection.onDocumentFormatting(this.onDocumentFormatting.bind(this));
        this.connection.onCodeAction(this.onCodeAction.bind(this));
        this.connection.languages.semanticTokens.on(this.onSemanticTokens.bind(this));
        this.connection.onDocumentSymbol(this.onDocumentSymbol.bind(this));
        this.connection.onFoldingRanges(this.onFoldingRanges.bind(this));

        this.documents.listen(this.connection);
    }

    private onInitialize(_params: InitializeParams): InitializeResult {
        logger.info('Initializing VASP Language Server (Modular)...');
        return {
            capabilities: {
                textDocumentSync: TextDocumentSyncKind.Incremental,
                hoverProvider: true,
                completionProvider: {
                    resolveProvider: true,
                    triggerCharacters: ['=']
                },
                documentFormattingProvider: true,
                codeActionProvider: true,
                documentSymbolProvider: true,
                semanticTokensProvider: {
                    legend: semanticTokensLegend,
                    full: true
                },
                foldingRangeProvider: true
            }
        };
    }

    private onInitialized(): void {
        logger.info('VASP Language Server Initialized.');
    }

    private async validateTextDocument(textDocument: TextDocument): Promise<void> {
        const uri = textDocument.uri;
        const diagnostics: Diagnostic[] = [];
        const fileName = uri.split('/').pop() || '';

        try {
            if (fileName.match(/POSCAR/i) || fileName.match(/CONTCAR/i)) {
                diagnostics.push(...validatePoscar(textDocument));
            } else if (fileName.match(/INCAR/i)) {
                const parsed = parseIncar(textDocument);
                diagnostics.push(...validateIncar(parsed));
            } else if (fileName.match(/POTCAR/i)) {
                diagnostics.push(...validatePotcar(textDocument));
            } else if (fileName.match(/KPOINTS/i)) {
                diagnostics.push(...validateKpoints(parseKpoints(textDocument)));
            }
        } catch (e) {
            logger.error(`Error validating ${fileName}: ${e}`);
        }

        this.connection.sendDiagnostics({ uri: textDocument.uri, diagnostics });
    }

    private onCompletion(_params: any): CompletionItem[] {
        const uri = _params.textDocument.uri;
        if (uri.match(/KPOINTS/i)) {
            return [
                {
                    label: 'Monkhorst-Pack',
                    kind: CompletionItemKind.Snippet,
                    insertText: 'Automatic Mesh\n0\nMonkhorst-Pack\n4 4 4\n0 0 0'
                },
                {
                    label: 'Gamma-Centered',
                    kind: CompletionItemKind.Snippet,
                    insertText: 'Automatic Mesh\n0\nGamma\n4 4 4\n0 0 0'
                }
            ];
        }
        return getIncarCompletions();
    }

    private onCompletionResolve(item: CompletionItem): CompletionItem {
        return resolveIncarCompletion(item);
    }

    private onHover({ textDocument, position }: any) {
        const document = this.documents.get(textDocument.uri);
        if (!document || !textDocument.uri.match(/INCAR/i)) return null;
        return getIncarHover(document, position);
    }

    private onDocumentFormatting(params: any) {
        const document = this.documents.get(params.textDocument.uri);
        if (document && document.uri.match(/INCAR/i)) {
            return formatIncar(document);
        }
        return [];
    }

    private onCodeAction(params: any) {
        return getIncarCodeActions(params.textDocument.uri, params.context.diagnostics);
    }

    private onSemanticTokens(params: any) {
        const document = this.documents.get(params.textDocument.uri);
        if (document && document.uri.match(/INCAR/i)) {
            return getIncarSemanticTokens(document);
        }
        return { data: [] };
    }

    private onDocumentSymbol(params: any) {
        const document = this.documents.get(params.textDocument.uri);
        if (!document) return null;
        if (document.uri.match(/INCAR/i)) return getIncarSymbols(document);
        if (document.uri.match(/POSCAR/i) || document.uri.match(/CONTCAR/i)) return getPoscarSymbols(document);
        return null;
    }

    private onFoldingRanges(params: any) {
        const document = this.documents.get(params.textDocument.uri);
        return document ? getFoldingRanges(document) : [];
    }

    /**
     * Starts the server.
     */
    public start(): void {
        this.connection.listen();
    }
}
