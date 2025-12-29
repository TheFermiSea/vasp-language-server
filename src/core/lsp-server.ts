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

// VASP Features
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

// CRYSTAL23 Features
import { parseCrystal } from '../features/crystal/parsing';
import { validateCrystal } from '../features/crystal/linting';
import { getCrystalCompletions, resolveCrystalCompletion } from '../features/crystal/completion';
import { getCrystalHover } from '../features/crystal/hover';
import { getCrystalSemanticTokens } from '../features/crystal/semantic-tokens';

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
        logger.info('Initializing DFT Language Server (VASP + CRYSTAL23)...');
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
        logger.info('DFT Language Server Initialized (VASP + CRYSTAL23).');
    }

    /**
     * Determine the file type from the document URI
     */
    private getFileType(uri: string): 'incar' | 'poscar' | 'kpoints' | 'potcar' | 'crystal' | 'unknown' {
        const fileName = uri.split('/').pop() || '';
        const upper = fileName.toUpperCase();

        // CRYSTAL23 .d12 files
        if (fileName.endsWith('.d12') || fileName.endsWith('.D12')) {
            return 'crystal';
        }

        // VASP files
        if (upper.includes('INCAR')) return 'incar';
        if (upper.includes('POSCAR') || upper.includes('CONTCAR')) return 'poscar';
        if (upper.includes('KPOINTS')) return 'kpoints';
        if (upper.includes('POTCAR')) return 'potcar';

        return 'unknown';
    }

    private async validateTextDocument(textDocument: TextDocument): Promise<void> {
        const uri = textDocument.uri;
        const diagnostics: Diagnostic[] = [];
        const fileType = this.getFileType(uri);

        try {
            let structure: VaspStructure | undefined;

            switch (fileType) {
                case 'crystal': {
                    const parsed = parseCrystal(textDocument);
                    structure = { type: 'crystal', data: parsed };
                    diagnostics.push(...validateCrystal(parsed));
                    break;
                }
                case 'poscar': {
                    const parsed = parsePoscar(textDocument);
                    structure = { type: 'poscar', data: parsed };
                    diagnostics.push(...validatePoscar(textDocument, parsed));
                    break;
                }
                case 'incar': {
                    const parsed = parseIncar(textDocument);
                    structure = { type: 'incar', data: parsed };
                    diagnostics.push(...validateIncar(parsed));
                    break;
                }
                case 'potcar': {
                    diagnostics.push(...(await validatePotcar(textDocument)));
                    break;
                }
                case 'kpoints': {
                    const parsed = parseKpoints(textDocument);
                    structure = { type: 'kpoints', data: parsed };
                    diagnostics.push(...validateKpoints(parsed));
                    break;
                }
            }

            if (structure) {
                this.cache.set(textDocument, structure);
            }
        } catch (e) {
            logger.error(`Error validating ${fileType} file: ${e}`);
        }

        this.connection.sendDiagnostics({ uri: textDocument.uri, diagnostics });
    }

    private onCompletion(_params: CompletionParams): CompletionItem[] {
        try {
            const document = this.documents.get(_params.textDocument.uri);
            if (!document) return [];

            const fileType = this.getFileType(document.uri);
            switch (fileType) {
                case 'crystal':
                    return getCrystalCompletions();
                case 'incar':
                    return getIncarCompletions();
                case 'kpoints':
                    return getKpointsCompletions(document.getText(), _params.position);
                default:
                    return [];
            }
        } catch (e) {
            logger.error(`Completion failed: ${e}`);
            return [];
        }
    }

    private onCompletionResolve(item: CompletionItem): CompletionItem {
        try {
            // Try CRYSTAL23 resolution first (check if it's a CRYSTAL keyword)
            const crystalResolved = resolveCrystalCompletion(item);
            if (crystalResolved.documentation) {
                return crystalResolved;
            }
            // Fall back to INCAR resolution
            return resolveIncarCompletion(item);
        } catch (e) {
            logger.error(`Error in onCompletionResolve: ${e}`);
            return item;
        }
    }

    private onHover(params: HoverParams) {
        try {
            const document = this.documents.get(params.textDocument.uri);
            if (!document) return null;

            const fileType = this.getFileType(document.uri);
            switch (fileType) {
                case 'crystal': {
                    const cached = this.cache.get(document);
                    if (cached?.type === 'crystal') {
                        return getCrystalHover(document, cached.data, params.position);
                    }
                    return getCrystalHover(document, parseCrystal(document), params.position);
                }
                case 'incar':
                    return getIncarHover(document, params.position);
                case 'poscar': {
                    const cached = this.cache.get(document);
                    if (cached?.type === 'poscar') {
                        return getPoscarHover(cached.data, params.position);
                    }
                    return getPoscarHover(parsePoscar(document), params.position);
                }
                case 'kpoints':
                    return getKpointsHover(document.getText(), params.position);
                default:
                    return null;
            }
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

            const fileType = this.getFileType(document.uri);
            switch (fileType) {
                case 'crystal': {
                    const cached = this.cache.get(document);
                    if (cached?.type === 'crystal') return getCrystalSemanticTokens(cached.data);
                    return getCrystalSemanticTokens(parseCrystal(document));
                }
                case 'incar': {
                    const cached = this.cache.get(document);
                    if (cached?.type === 'incar') return getIncarSemanticTokens(document, cached.data);
                    return getIncarSemanticTokens(document, parseIncar(document));
                }
                case 'poscar': {
                    const cached = this.cache.get(document);
                    if (cached?.type === 'poscar') return getPoscarSemanticTokens(cached.data);
                    return getPoscarSemanticTokens(parsePoscar(document));
                }
                case 'kpoints':
                    return getKpointsSemanticTokens(document.getText());
                default:
                    return { data: [] };
            }
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
