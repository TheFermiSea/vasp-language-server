import { Hover, Position } from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { parseIncar } from './parsing';
import { VASP_TAGS } from '../../data/vasp-tags';

/**
 * Provide hover documentation for INCAR tags and values.
 *
 * @param document - LSP text document for an INCAR file.
 * @param position - Hover position.
 * @returns Hover payload with tag documentation, or null if none.
 */
export function getIncarHover(document: TextDocument, position: Position): Hover | null {
    const parsed = parseIncar(document);
    const line = position.line;
    const char = position.character;

    const token = parsed.allTokens.find((t) => {
        return t.range.start.line === line && t.range.start.character <= char && t.range.end.character >= char;
    });

    if (token) {
        if (token.type === 'value' || token.type === 'tag') {
            const tagName = token.text.toUpperCase();
            const def = VASP_TAGS[tagName];

            if (def) {
                const markup = `**${tagName}**\n\n${def.description}\n\n*Type*: \`${def.type}\`\n*Default*: ${def.default || 'N/A'}`;
                return {
                    contents: {
                        kind: 'markdown',
                        value: markup
                    }
                };
            }
        }
    }
    return null;
}
