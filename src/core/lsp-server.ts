import {
    createConnection,
    TextDocuments,
    ProposedFeatures,
    InitializeParams,
    TextDocumentSyncKind,
    InitializeResult,
    Diagnostic,
    CompletionItem,
    _Connection,
    CompletionParams,
    HoverParams,
    DocumentFormattingParams,
    CodeActionParams,
    SemanticTokensParams,
    DocumentSymbolParams,
    FoldingRangeParams
} from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import * as path from 'path';

// Features
import { validatePoscar } from '../features/poscar/linting';
import { validateIncar } from '../features/incar/linting';
import { validatePotcar } from '../features/potcar/linting';
import { parseKpoints } from '../features/kpoints/parsing';
import { validateKpoints } from '../features/kpoints/linting';
import { parseIncar } from '../features/incar/parsing';
import { parsePoscar } from '../features/poscar/parsing';
import { formatIncar } from '../features/incar/formatting';
import { getIncarCodeActions } from '../features/incar/code-actions';
import { getIncarSymbols } from '../features/incar/symbols';
import { getPoscarSymbols } from '../features/poscar/symbols';
import { getIncarHover } from '../features/incar/hover';
import { getIncarCompletions, resolveIncarCompletion } from '../features/incar/completion';
import { getIncarSemanticTokens } from '../features/incar/semantic-tokens';
import { getPoscarHover } from '../features/poscar/hover';
import { getPoscarSemanticTokens } from '../features/poscar/semantic-tokens';
import { getKpointsCompletions } from '../features/kpoints/completion';
import { getKpointsHover } from '../features/kpoints/hover';
import { getKpointsSemanticTokens } from '../features/kpoints/semantic-tokens';
import { semanticTokensLegend } from '../features/semantic-tokens-legend';
import { getFoldingRanges } from '../features/poscar/folding';
import { getIncarFoldingRanges } from '../features/incar/folding';

// Core & Utils
import { logger } from '../utils/logger';
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
            const document = this.documents.get(_params.textDocument.uri);
            if (!document) return [];

            const basename = path.basename(document.uri).toUpperCase();
            if (basename.includes('INCAR')) {
                return getIncarCompletions();
            } else if (basename.includes('KPOINTS')) {
                return getKpointsCompletions(document.getText(), _params.position);
            }
            return [];
        } catch (e) {
            logger.error(`Completion failed: ${e}`);
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

            const basename = path.basename(document.uri).toUpperCase();
            if (basename.includes('INCAR')) {
                return getIncarHover(document, params.position);
            } else if (basename.includes('POSCAR') || basename.includes('CONTCAR')) {
                const cached = this.cache.get(document);
                if (cached?.type === 'poscar') {
                    return getPoscarHover(cached.data, params.position);
                } else {
                    return getPoscarHover(parsePoscar(document), params.position);
                }
            } else if (basename.includes('KPOINTS')) {
                return getKpointsHover(document.getText(), params.position);
            }
            return null;
        } catch (e) {
            logger.error(`Hover failed: ${e}`);
            return null;
        }
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

            const basename = path.basename(document.uri).toUpperCase();
            if (basename.includes('INCAR')) {
                const cached = this.cache.get(document);
                if (cached?.type === 'incar') return getIncarSemanticTokens(document, cached.data);
                return getIncarSemanticTokens(document, parseIncar(document));
            } else if (basename.includes('POSCAR') || basename.includes('CONTCAR')) {
                const cached = this.cache.get(document);
                if (cached?.type === 'poscar') return getPoscarSemanticTokens(cached.data);
                return getPoscarSemanticTokens(parsePoscar(document));
            } else if (basename.includes('KPOINTS')) {
                return getKpointsSemanticTokens(document.getText());
            }
            return { data: [] };
        } catch (e) {
            logger.error(`SemanticTokens failed: ${e}`);
            return { data: [] };
        }
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
            if (document.uri.match(/POSCAR/i) || document.uri.match(/CONTCAR/i))
                return getPoscarSymbols(document, parsePoscar(document).lines);
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
            if (document.uri.match(/INCAR/i)) {
                return getIncarFoldingRanges(document);
            }
            const lines = document.getText().split(/\r?\n/);
            return getFoldingRanges(lines, cached?.type === 'poscar' ? cached.data.lines : undefined);
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
