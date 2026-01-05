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
 * Configuration options for the LSP server.
 */
export interface LspServerOptions {
    /** Maximum number of documents to cache (default: 50) */
    cacheSize?: number;
}

/**
 * Encapsulates the VASP Language Server logic.
 */
export class LspServer {
    private connection: _Connection;
    private documents: TextDocuments<TextDocument>;
    private cache: DocumentCache;

    constructor(options: LspServerOptions = {}) {
        this.connection = createConnection(ProposedFeatures.all);
        this.documents = new TextDocuments(TextDocument);
        this.cache = new DocumentCache(options.cacheSize);

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
     * Determine the file type from the document URI.
     *
     * Handles various edge cases:
     * - Cross-platform path separators (/ and \)
     * - URI-encoded characters (e.g., %20 for spaces)
     * - Case-insensitive matching for VASP files
     * - File extensions (.d12 for CRYSTAL23)
     * - Common naming patterns (INCAR, POSCAR, CONTCAR, KPOINTS, POTCAR)
     * - Prefixed/suffixed filenames (e.g., my_INCAR, INCAR.bak, INCAR_relaxation)
     */
    private getFileType(uri: string): 'incar' | 'poscar' | 'kpoints' | 'potcar' | 'crystal' | 'unknown' {
        // Decode URI-encoded characters (e.g., %20 -> space)
        let decodedUri: string;
        try {
            decodedUri = decodeURIComponent(uri);
        } catch {
            // If decoding fails, use original URI
            decodedUri = uri;
        }

        // Extract filename - handle both forward and back slashes for cross-platform support
        // Also handle 'file://' protocol prefix
        const pathWithoutProtocol = decodedUri.replace(/^file:\/\//, '');
        const pathParts = pathWithoutProtocol.split(/[/\\]/);
        const fileName = pathParts[pathParts.length - 1] || '';

        // Handle empty filename edge case
        if (!fileName.trim()) {
            return 'unknown';
        }

        // Normalize to uppercase for case-insensitive comparison
        const upperFileName = fileName.toUpperCase();

        // CRYSTAL23 .d12 files - check extension case-insensitively
        if (upperFileName.endsWith('.D12')) {
            return 'crystal';
        }

        // VASP files - use word boundary-aware matching to avoid false positives
        // This handles: INCAR, my_INCAR, INCAR.bak, INCAR_old, but not MINCAR

        // Check for INCAR (but not things like MINCAR - require word boundary or start)
        if (this.matchesVaspFilePattern(upperFileName, 'INCAR')) {
            return 'incar';
        }

        // Check for POSCAR or CONTCAR
        if (
            this.matchesVaspFilePattern(upperFileName, 'POSCAR') ||
            this.matchesVaspFilePattern(upperFileName, 'CONTCAR')
        ) {
            return 'poscar';
        }

        // Check for KPOINTS
        if (this.matchesVaspFilePattern(upperFileName, 'KPOINTS')) {
            return 'kpoints';
        }

        // Check for POTCAR
        if (this.matchesVaspFilePattern(upperFileName, 'POTCAR')) {
            return 'potcar';
        }

        return 'unknown';
    }

    /**
     * Check if a filename matches a VASP file pattern.
     *
     * Matches patterns like:
     * - Exact match: INCAR
     * - With prefix: my_INCAR, test-INCAR
     * - With suffix: INCAR.bak, INCAR_old, INCAR.1
     * - Combined: my_INCAR.bak
     *
     * Does NOT match if the pattern is embedded within a larger word:
     * - MINCAR (M is a letter, not a valid prefix separator)
     *
     * @param upperFileName - The uppercase filename to check
     * @param pattern - The VASP file pattern to match (e.g., 'INCAR')
     */
    private matchesVaspFilePattern(upperFileName: string, pattern: string): boolean {
        const index = upperFileName.indexOf(pattern);
        if (index === -1) {
            return false;
        }

        // Check character before the pattern (if any)
        // Valid: start of string, or non-alphanumeric character (like _, -, .)
        if (index > 0) {
            const charBefore = upperFileName[index - 1];
            if (/[A-Z0-9]/.test(charBefore)) {
                return false;
            }
        }

        // Check character after the pattern (if any)
        // Valid: end of string, or non-alphanumeric character (like _, -, ., or extension)
        const afterIndex = index + pattern.length;
        if (afterIndex < upperFileName.length) {
            const charAfter = upperFileName[afterIndex];
            if (/[A-Z0-9]/.test(charAfter)) {
                return false;
            }
        }

        return true;
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
        } catch (error) {
            logger.error(`Validation failed for ${fileType} file`, error);
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
        } catch (error) {
            logger.error('Completion failed', error);
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
        } catch (error) {
            logger.error('Completion resolve failed', error);
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
        } catch (error) {
            logger.error('Hover failed', error);
            return null;
        }
    }

    private onDocumentFormatting(params: DocumentFormattingParams) {
        try {
            const document = this.documents.get(params.textDocument.uri);
            if (!document) return [];

            const fileType = this.getFileType(document.uri);
            if (fileType === 'incar') {
                return formatIncar(document);
            }
        } catch (error) {
            logger.error('Document formatting failed', error);
        }
        return [];
    }

    private onCodeAction(params: CodeActionParams) {
        try {
            const fileType = this.getFileType(params.textDocument.uri);
            if (fileType === 'incar') {
                return getIncarCodeActions(params.textDocument.uri, params.context.diagnostics);
            }
            return [];
        } catch (error) {
            logger.error('Code action failed', error);
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
        } catch (error) {
            logger.error('Semantic tokens failed', error);
            return { data: [] };
        }
    }

    private onDocumentSymbol(params: DocumentSymbolParams) {
        try {
            const document = this.documents.get(params.textDocument.uri);
            if (!document) return null;

            const fileType = this.getFileType(document.uri);
            const cached = this.cache.get(document);

            switch (fileType) {
                case 'incar': {
                    if (cached?.type === 'incar') return getIncarSymbols(document, cached.data);
                    return getIncarSymbols(document, parseIncar(document));
                }
                case 'poscar': {
                    if (cached?.type === 'poscar') return getPoscarSymbols(document, cached.data.lines);
                    return getPoscarSymbols(document, parsePoscar(document).lines);
                }
                default:
                    return null;
            }
        } catch (error) {
            logger.error('Document symbol failed', error);
        }
        return null;
    }

    private onFoldingRanges(params: FoldingRangeParams) {
        try {
            const document = this.documents.get(params.textDocument.uri);
            if (!document) return [];

            const fileType = this.getFileType(document.uri);
            const cached = this.cache.get(document);

            switch (fileType) {
                case 'incar':
                    return getIncarFoldingRanges(document);
                case 'poscar': {
                    const lines = document.getText().split(/\r?\n/);
                    const poscarLines = cached?.type === 'poscar' ? cached.data.lines : undefined;
                    return getFoldingRanges(lines, poscarLines);
                }
                default:
                    return [];
            }
        } catch (error) {
            logger.error('Folding ranges failed', error);
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
