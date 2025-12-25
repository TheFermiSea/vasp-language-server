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
    _Connection,
    CompletionParams, // ADDED
    HoverParams, // ADDED
    DocumentFormattingParams, // ADDED
    CodeActionParams, // ADDED
    SemanticTokensParams, // ADDED
    DocumentSymbolParams, // ADDED
    FoldingRangeParams // ADDED
} from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { validatePoscar } from './poscar-linting';
import { validateIncar } from './incar-linting';
import { validatePotcar } from './potcar-linting';
import { parseKpoints } from './kpoints-parsing';
import { validateKpoints } from './kpoints-linting';
import { parseIncar, IncarDocument } from './incar-parsing';
import { parsePoscar, PoscarDocument } from './poscar-parsing';
import { formatIncar } from './incar-formatting';
import { logger } from './logger';
import { getIncarCodeActions } from './code-actions';
import { getIncarSymbols, getPoscarSymbols } from './document-symbols';
import { getIncarHover } from './features/incar/hover';
import { getIncarCompletions, resolveIncarCompletion } from './features/incar/completion';
import { getIncarSemanticTokens } from './features/incar/semantic-tokens';
import { semanticTokensLegend } from './features/semantic-tokens-legend';
import { getFoldingRanges } from './features/folding';
import { DocumentCache, VaspStructure } from './document-cache';

/**
 * Encapsulates the VASP Language Server logic.
 */
export class LspServer {
    private connection: _Connection;
    private documents: TextDocuments<TextDocument>;
    private cache: DocumentCache;

    constructor() {
        this.connection = createConnection(ProposedFeatures.all);
        this.documents = new TextDocuments(TextDocument);
        this.cache = new DocumentCache();

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

        this.documents.onDidClose((e) => {
            this.cache.delete(e.document.uri);
        });

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
            let structure: VaspStructure | undefined;

            if (fileName.match(/POSCAR/i) || fileName.match(/CONTCAR/i)) {
                const parsed = parsePoscar(textDocument);
                structure = { type: 'poscar', data: parsed };
                diagnostics.push(...validatePoscar(textDocument));
            } else if (fileName.match(/INCAR/i)) {
                const parsed = parseIncar(textDocument);
                structure = { type: 'incar', data: parsed };
                diagnostics.push(...validateIncar(parsed));
            } else if (fileName.match(/POTCAR/i)) {
                diagnostics.push(...(await validatePotcar(textDocument)));
            } else if (fileName.match(/KPOINTS/i)) {
                const parsed = parseKpoints(textDocument);
                structure = { type: 'kpoints', data: parsed };
                diagnostics.push(...validateKpoints(parsed));
            }

            if (structure) {
                this.cache.set(textDocument, structure);
            }
        } catch (e) {
            logger.error(`Error validating ${fileName}: ${e} `);
        }

        this.connection.sendDiagnostics({ uri: textDocument.uri, diagnostics });
    }

    private onCompletion(_params: CompletionParams): CompletionItem[] {
        try {
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
        } catch (e) {
            logger.error(`Error in onCompletion: ${e} `);
            return [];
        }
    }

    private onCompletionResolve(item: CompletionItem): CompletionItem {
        try {
            return resolveIncarCompletion(item);
        } catch (e) {
            logger.error(`Error in onCompletionResolve: ${e} `);
            return item;
        }
    }

    private onHover(params: HoverParams) {
        try {
            const document = this.documents.get(params.textDocument.uri);
            if (!document) return null;

            const cached = this.cache.get(document);
            if (cached?.type === 'incar') {
                return getIncarHover(document, params.position);
            }
        } catch (e) {
            logger.error(`Error in onHover: ${e} `);
        }
        return null;
    }

    private onDocumentFormatting(params: DocumentFormattingParams) {
        try {
            const document = this.documents.get(params.textDocument.uri);
            if (document && document.uri.match(/INCAR/i)) {
                return formatIncar(document);
            }
        } catch (e) {
            logger.error(`Error in onDocumentFormatting: ${e} `);
        }
        return [];
    }

    private onCodeAction(params: CodeActionParams) {
        try {
            return getIncarCodeActions(params.textDocument.uri, params.context.diagnostics);
        } catch (e) {
            logger.error(`Error in onCodeAction: ${e} `);
            return [];
        }
    }

    private onSemanticTokens(params: SemanticTokensParams) {
        try {
            const document = this.documents.get(params.textDocument.uri);
            if (!document) return { data: [] };

            const cached = this.cache.get(document);
            if (cached?.type === 'incar') {
                return getIncarSemanticTokens(document, cached.data);
            }
        } catch (e) {
            logger.error(`Error in onSemanticTokens: ${e} `);
        }
        return { data: [] };
    }

    private onDocumentSymbol(params: DocumentSymbolParams) {
        try {
            const document = this.documents.get(params.textDocument.uri);
            if (!document) return null;

            const cached = this.cache.get(document);
            if (cached?.type === 'incar') return getIncarSymbols(document, cached.data);
            if (cached?.type === 'poscar') return getPoscarSymbols(document, cached.data.lines);

            // Fallback or specific check for files not cached yet
            if (document.uri.match(/INCAR/i)) return getIncarSymbols(document, parseIncar(document));
            if (document.uri.match(/POSCAR/i) || document.uri.match(/CONTCAR/i)) return getPoscarSymbols(document, parsePoscar(document).lines);
        } catch (e) {
            logger.error(`Error in onDocumentSymbol: ${e} `);
        }
        return null;
    }

    private onFoldingRanges(params: FoldingRangeParams) {
        try {
            const document = this.documents.get(params.textDocument.uri);
            if (!document) return [];

            const cached = this.cache.get(document);
            return getFoldingRanges(document, cached?.type === 'poscar' ? cached.data.lines : undefined);
        } catch (e) {
            logger.error(`Error in onFoldingRanges: ${e} `);
            return [];
        }
    }

    /**
     * Starts the server.
     */
    public start(): void {
        this.connection.listen();
    }
}
