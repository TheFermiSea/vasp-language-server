/**
 * CRYSTAL23 semantic token highlighting
 */

import { SemanticTokensBuilder } from 'vscode-languageserver';
import { CrystalDocument, CrystalToken } from './parsing';
import { TokenType, TokenModifier } from '../semantic-tokens-legend';

/**
 * Generate semantic tokens for a CRYSTAL23 document.
 *
 * @param document - Parsed CRYSTAL document.
 * @returns Semantic tokens payload for the LSP response.
 */
export function getCrystalSemanticTokens(document: CrystalDocument): { data: number[] } {
    const builder = new SemanticTokensBuilder();

    // Sort tokens by position
    const sortedTokens = [...document.allTokens].sort((a, b) => {
        if (a.range.start.line !== b.range.start.line) {
            return a.range.start.line - b.range.start.line;
        }
        return a.range.start.character - b.range.start.character;
    });

    for (const token of sortedTokens) {
        const { tokenType, tokenModifiers } = getSemanticTokenInfo(token);
        if (tokenType !== undefined) {
            const length = token.range.end.character - token.range.start.character;
            builder.push(token.range.start.line, token.range.start.character, length, tokenType, tokenModifiers);
        }
    }

    return builder.build();
}

function getSemanticTokenInfo(token: CrystalToken): { tokenType?: number; tokenModifiers: number } {
    let tokenType: number | undefined;
    let tokenModifiers = 0;

    switch (token.type) {
        case 'title':
            tokenType = TokenType.comment;
            break;

        case 'geometry-type':
            tokenType = TokenType.class;
            tokenModifiers = 1 << TokenModifier.declaration;
            break;

        case 'keyword':
            tokenType = TokenType.keyword;
            break;

        case 'block-start':
            tokenType = TokenType.keyword;
            tokenModifiers = 1 << TokenModifier.declaration;
            break;

        case 'block-end':
        case 'end':
            tokenType = TokenType.keyword;
            break;

        case 'space-group':
        case 'atom-count':
            tokenType = TokenType.number;
            tokenModifiers = 1 << TokenModifier.declaration;
            break;

        case 'lattice-param':
        case 'number':
            tokenType = TokenType.number;
            break;

        case 'atom-line':
            tokenType = TokenType.variable;
            break;

        case 'basis-header':
            tokenType = TokenType.type;
            break;

        case 'basis-data':
            tokenType = TokenType.number;
            break;

        case 'comment':
            tokenType = TokenType.comment;
            break;

        case 'unknown':
            // Don't highlight unknown tokens
            break;
    }

    return { tokenType, tokenModifiers };
}
