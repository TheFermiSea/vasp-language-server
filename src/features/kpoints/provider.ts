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
import { parseKpoints } from './parsing';
import { validateKpoints } from './linting';
import { getKpointsCompletions } from './completion';
import { getKpointsHover } from './hover';
import { getKpointsSemanticTokens } from './semantic-tokens';

export class KpointsFeatureProvider extends BaseFeatureProvider {
    parse(document: TextDocument): VaspStructure {
        return { type: 'kpoints', data: parseKpoints(document) };
    }

    async validate(_document: TextDocument, parsed?: VaspStructure): Promise<Diagnostic[]> {
        const kpoints = parsed?.type === 'kpoints' ? parsed.data : undefined;
        return validateKpoints(kpoints ?? parseKpoints(_document));
    }

    complete(document: TextDocument, params: CompletionParams): CompletionItem[] {
        return getKpointsCompletions(document.getText(), params.position);
    }

    hover(document: TextDocument, params: HoverParams): Hover | null {
        return getKpointsHover(document.getText(), params.position);
    }

    getSemanticTokens(document: TextDocument, _params: SemanticTokensParams): SemanticTokens {
        return getKpointsSemanticTokens(document.getText());
    }
}
