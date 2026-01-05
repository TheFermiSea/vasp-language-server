import {
    Hover,
    HoverParams,
    SemanticTokens,
    SemanticTokensParams,
    DocumentSymbol,
    DocumentSymbolParams,
    FoldingRange,
    FoldingRangeParams,
    Diagnostic
} from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { BaseFeatureProvider } from '../../core/feature-provider';
import { VaspStructure } from '../../core/document-cache';
import { parsePoscar } from './parsing';
import { validatePoscar } from './linting';
import { getPoscarHover } from './hover';
import { getPoscarSemanticTokens } from './semantic-tokens';
import { getPoscarSymbols } from './symbols';
import { getFoldingRanges } from './folding';

export class PoscarFeatureProvider extends BaseFeatureProvider {
    parse(document: TextDocument): VaspStructure {
        return { type: 'poscar', data: parsePoscar(document) };
    }

    async validate(document: TextDocument, parsed?: VaspStructure): Promise<Diagnostic[]> {
        const poscar = parsed?.type === 'poscar' ? parsed.data : undefined;
        return validatePoscar(document, poscar ?? parsePoscar(document));
    }

    hover(_document: TextDocument, params: HoverParams, parsed?: VaspStructure): Hover | null {
        const poscar = parsed?.type === 'poscar' ? parsed.data : undefined;
        return getPoscarHover(poscar ?? parsePoscar(_document), params.position);
    }

    getSemanticTokens(_document: TextDocument, _params: SemanticTokensParams, parsed?: VaspStructure): SemanticTokens {
        const poscar = parsed?.type === 'poscar' ? parsed.data : undefined;
        return getPoscarSemanticTokens(poscar ?? parsePoscar(_document));
    }

    getSymbols(document: TextDocument, _params: DocumentSymbolParams, parsed?: VaspStructure): DocumentSymbol[] | null {
        const poscar = parsed?.type === 'poscar' ? parsed.data : undefined;
        return getPoscarSymbols(document, (poscar ?? parsePoscar(document)).lines);
    }

    getFoldingRanges(document: TextDocument, _params: FoldingRangeParams, parsed?: VaspStructure): FoldingRange[] {
        const lines = document.getText().split(/\r?\n/);
        const poscarLines = parsed?.type === 'poscar' ? parsed.data.lines : undefined;
        return getFoldingRanges(lines, poscarLines);
    }
}
