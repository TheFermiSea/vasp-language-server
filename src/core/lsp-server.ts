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

import { resolveIncarCompletion } from '../features/incar/completion';
import { resolveCrystalCompletion } from '../features/crystal/completion';
import { semanticTokensLegend } from '../features/semantic-tokens-legend';

// Core & Utils
import { logger } from '../utils/logger';
import { DocumentCache } from './document-cache';
import { FeatureProvider, FileType } from './feature-provider';
import { IncarFeatureProvider } from '../features/incar/provider';
import { PoscarFeatureProvider } from '../features/poscar/provider';
import { KpointsFeatureProvider } from '../features/kpoints/provider';
import { PotcarFeatureProvider } from '../features/potcar/provider';
import { CrystalFeatureProvider } from '../features/crystal/provider';

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
    private providers: Map<FileType, FeatureProvider>;

    /**
     * Create a new language server instance with optional configuration.
     *
     * @param options - Server configuration options.
     */
    constructor(options: LspServerOptions = {}) {
        this.connection = createConnection(ProposedFeatures.all);
        this.documents = new TextDocuments(TextDocument);
        this.cache = new DocumentCache(options.cacheSize);
        this.providers = new Map<FileType, FeatureProvider>([
            ['incar', new IncarFeatureProvider()],
            ['poscar', new PoscarFeatureProvider()],
            ['kpoints', new KpointsFeatureProvider()],
            ['potcar', new PotcarFeatureProvider()],
            ['crystal', new CrystalFeatureProvider()]
        ]);

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
    private getFileType(uri: string): FileType {
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
        const provider = this.providers.get(fileType);

        try {
            if (provider) {
                const structure = provider.parse ? provider.parse(textDocument) : undefined;
                if (structure) {
                    this.cache.set(textDocument, structure);
                }
                diagnostics.push(...(await provider.validate(textDocument, structure)));
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
            const provider = this.providers.get(fileType);
            return provider ? provider.complete(document, _params) : [];
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
            const provider = this.providers.get(fileType);
            const cached = this.cache.get(document);
            return provider ? provider.hover(document, params, cached) : null;
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
            const provider = this.providers.get(fileType);
            return provider ? provider.format(document, params) : [];
        } catch (error) {
            logger.error('Document formatting failed', error);
        }
        return [];
    }

    private onCodeAction(params: CodeActionParams) {
        try {
            const document = this.documents.get(params.textDocument.uri);
            if (!document) return [];
            const fileType = this.getFileType(params.textDocument.uri);
            const provider = this.providers.get(fileType);
            return provider ? provider.getCodeActions(document, params) : [];
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
            const provider = this.providers.get(fileType);
            const cached = this.cache.get(document);
            return provider ? provider.getSemanticTokens(document, params, cached) : { data: [] };
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
            const provider = this.providers.get(fileType);
            const cached = this.cache.get(document);
            return provider ? provider.getSymbols(document, params, cached) : null;
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
            const provider = this.providers.get(fileType);
            const cached = this.cache.get(document);
            return provider ? provider.getFoldingRanges(document, params, cached) : [];
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
