import {
    CodeAction,
    CodeActionParams,
    CompletionItem,
    CompletionParams,
    Diagnostic,
    DocumentFormattingParams,
    DocumentSymbol,
    DocumentSymbolParams,
    FoldingRange,
    FoldingRangeParams,
    Hover,
    HoverParams,
    SemanticTokens,
    SemanticTokensParams,
    TextEdit
} from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { VaspStructure } from './document-cache';

export type FileType = 'incar' | 'poscar' | 'kpoints' | 'potcar' | 'crystal' | 'unknown';

export interface FeatureProvider {
    parse?(document: TextDocument): VaspStructure | undefined;
    validate(document: TextDocument, parsed?: VaspStructure): Promise<Diagnostic[]>;
    complete(document: TextDocument, params: CompletionParams): CompletionItem[];
    hover(document: TextDocument, params: HoverParams, parsed?: VaspStructure): Hover | null;
    format(document: TextDocument, params: DocumentFormattingParams): TextEdit[];
    getSemanticTokens(document: TextDocument, params: SemanticTokensParams, parsed?: VaspStructure): SemanticTokens;
    getSymbols(document: TextDocument, params: DocumentSymbolParams, parsed?: VaspStructure): DocumentSymbol[] | null;
    getFoldingRanges(document: TextDocument, params: FoldingRangeParams, parsed?: VaspStructure): FoldingRange[];
    getCodeActions(document: TextDocument, params: CodeActionParams): CodeAction[];
}

export abstract class BaseFeatureProvider implements FeatureProvider {
    parse(_document: TextDocument): VaspStructure | undefined {
        return undefined;
    }

    async validate(_document: TextDocument, _parsed?: VaspStructure): Promise<Diagnostic[]> {
        return [];
    }

    complete(_document: TextDocument, _params: CompletionParams): CompletionItem[] {
        return [];
    }

    hover(_document: TextDocument, _params: HoverParams, _parsed?: VaspStructure): Hover | null {
        return null;
    }

    format(_document: TextDocument, _params: DocumentFormattingParams): TextEdit[] {
        return [];
    }

    getSemanticTokens(_document: TextDocument, _params: SemanticTokensParams, _parsed?: VaspStructure): SemanticTokens {
        return { data: [] };
    }

    getSymbols(
        _document: TextDocument,
        _params: DocumentSymbolParams,
        _parsed?: VaspStructure
    ): DocumentSymbol[] | null {
        return null;
    }

    getFoldingRanges(_document: TextDocument, _params: FoldingRangeParams, _parsed?: VaspStructure): FoldingRange[] {
        return [];
    }

    getCodeActions(_document: TextDocument, _params: CodeActionParams): CodeAction[] {
        return [];
    }
}
