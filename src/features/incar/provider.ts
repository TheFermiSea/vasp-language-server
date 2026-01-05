import {
    CompletionItem,
    CompletionParams,
    DocumentFormattingParams,
    Hover,
    HoverParams,
    SemanticTokens,
    SemanticTokensParams,
    DocumentSymbol,
    DocumentSymbolParams,
    FoldingRange,
    FoldingRangeParams,
    CodeAction,
    CodeActionParams,
    TextEdit,
    Diagnostic
} from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { BaseFeatureProvider } from '../../core/feature-provider';
import { VaspStructure } from '../../core/document-cache';
import { parseIncar } from './parsing';
import { validateIncar } from './linting';
import { getIncarCompletions } from './completion';
import { getIncarHover } from './hover';
import { getIncarSemanticTokens } from './semantic-tokens';
import { getIncarSymbols } from './symbols';
import { getIncarFoldingRanges } from './folding';
import { formatIncar } from './formatting';
import { getIncarCodeActions } from './code-actions';

export class IncarFeatureProvider extends BaseFeatureProvider {
    parse(document: TextDocument): VaspStructure {
        return { type: 'incar', data: parseIncar(document) };
    }

    async validate(_document: TextDocument, parsed?: VaspStructure): Promise<Diagnostic[]> {
        const incar = parsed?.type === 'incar' ? parsed.data : undefined;
        return validateIncar(incar ?? parseIncar(_document));
    }

    complete(_document: TextDocument, _params: CompletionParams): CompletionItem[] {
        return getIncarCompletions();
    }

    hover(document: TextDocument, params: HoverParams): Hover | null {
        return getIncarHover(document, params.position);
    }

    format(document: TextDocument, _params: DocumentFormattingParams): TextEdit[] {
        return formatIncar(document);
    }

    getSemanticTokens(document: TextDocument, _params: SemanticTokensParams, parsed?: VaspStructure): SemanticTokens {
        const incar = parsed?.type === 'incar' ? parsed.data : undefined;
        return getIncarSemanticTokens(document, incar ?? parseIncar(document));
    }

    getSymbols(document: TextDocument, _params: DocumentSymbolParams, parsed?: VaspStructure): DocumentSymbol[] | null {
        const incar = parsed?.type === 'incar' ? parsed.data : undefined;
        return getIncarSymbols(document, incar ?? parseIncar(document));
    }

    getFoldingRanges(document: TextDocument, _params: FoldingRangeParams): FoldingRange[] {
        return getIncarFoldingRanges(document);
    }

    getCodeActions(_document: TextDocument, params: CodeActionParams): CodeAction[] {
        return getIncarCodeActions(params.textDocument.uri, params.context.diagnostics);
    }
}
