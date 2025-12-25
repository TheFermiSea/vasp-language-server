import { SemanticTokensBuilder, SemanticTokens } from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { parseIncar, IncarDocument } from '../../incar-parsing';
import { TokenType, TokenModifier } from '../semantic-tokens-legend';
import { isNumber } from '../../util';

export function getIncarSemanticTokens(document: TextDocument, parsed: IncarDocument): SemanticTokens {
    const builder = new SemanticTokensBuilder();

    // Collect all items to simplify sorting
    interface HighlightItem {
        line: number;
        char: number;
        length: number;
        type: number;
        modifiers: number;
    }

    const items: HighlightItem[] = [];

    // 1. Tags and Values from Statements
    for (const stmt of parsed.statements) {
        // Tag -> 'property' + 'declaration'
        const tag = stmt.tag;
        items.push({
            line: tag.range.start.line,
            char: tag.range.start.character,
            length: tag.text.length,
            type: TokenType.property,
            modifiers: 1 << TokenModifier.declaration
        });

        // Values
        for (const val of stmt.values) {
            let type = TokenType.string;
            const text = val.text;

            if (isNumber(text)) {
                type = TokenType.number;
            } else {
                const lower = text.toLowerCase();
                if (
                    lower === '.true.' ||
                    lower === '.false.' ||
                    lower === 't' ||
                    lower === 'f' ||
                    lower === 'true' ||
                    lower === 'false'
                ) {
                    type = TokenType.keyword;
                }
            }

            items.push({
                line: val.range.start.line,
                char: val.range.start.character,
                length: val.text.length,
                type: type,
                modifiers: 0
            });
        }
    }

    // 2. Comments from allTokens
    for (const t of parsed.allTokens) {
        if (t.type === 'comment') {
            items.push({
                line: t.range.start.line,
                char: t.range.start.character,
                length: t.text.length,
                type: TokenType.comment,
                modifiers: 0
            });
        }
    }

    // 3. Sort by position
    items.sort((a, b) => {
        if (a.line !== b.line) {
            return a.line - b.line;
        }
        return a.char - b.char;
    });

    // 4. Push to builder
    for (const item of items) {
        builder.push(item.line, item.char, item.length, item.type, item.modifiers);
    }

    return builder.build();
}
