import { Hover, Position } from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { parseIncar } from '../../incar-parsing';
import { VASP_TAGS } from '../../data/vasp-tags';

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
