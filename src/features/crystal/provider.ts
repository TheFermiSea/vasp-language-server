import {
    CompletionItem,
    CompletionParams,
    Hover,
    HoverParams,
    SemanticTokens,
    SemanticTokensParams,
    Diagnostic
} from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { BaseFeatureProvider } from '../../core/feature-provider';
import { VaspStructure } from '../../core/document-cache';
import { parseCrystal } from './parsing';
import { validateCrystal } from './linting';
import { getCrystalCompletions } from './completion';
import { getCrystalHover } from './hover';
import { getCrystalSemanticTokens } from './semantic-tokens';

export class CrystalFeatureProvider extends BaseFeatureProvider {
    parse(document: TextDocument): VaspStructure {
        return { type: 'crystal', data: parseCrystal(document) };
    }

    async validate(_document: TextDocument, parsed?: VaspStructure): Promise<Diagnostic[]> {
        const crystal = parsed?.type === 'crystal' ? parsed.data : undefined;
        return validateCrystal(crystal ?? parseCrystal(_document));
    }

    complete(_document: TextDocument, _params: CompletionParams): CompletionItem[] {
        return getCrystalCompletions();
    }

    hover(document: TextDocument, params: HoverParams, parsed?: VaspStructure): Hover | null {
        const crystal = parsed?.type === 'crystal' ? parsed.data : undefined;
        return getCrystalHover(document, crystal ?? parseCrystal(document), params.position);
    }

    getSemanticTokens(_document: TextDocument, _params: SemanticTokensParams, parsed?: VaspStructure): SemanticTokens {
        const crystal = parsed?.type === 'crystal' ? parsed.data : undefined;
        return getCrystalSemanticTokens(crystal ?? parseCrystal(_document));
    }
}
