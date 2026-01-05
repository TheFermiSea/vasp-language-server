import { TextDocument } from 'vscode-languageserver-textdocument';
import { Position, Range, TextEdit } from 'vscode-languageserver-types';
import { parsePoscar } from './parsing';

/**
 * Format a POSCAR document by normalizing whitespace.
 *
 * @param document - LSP text document for a POSCAR/CONTCAR file.
 * @returns Text edits that normalize whitespace for each line.
 */
export function formatPoscar(document: TextDocument): TextEdit[] {
    const parsed = parsePoscar(document);
    const edits: TextEdit[] = [];

    for (const line of parsed.lines) {
        const original = line.line.text;
        if (line.tokens.length === 0) {
            continue;
        }

        let next = original;
        if (line.type === 'comment') {
            next = original.replace(/\s+$/, '');
        } else {
            next = line.tokens.map((token) => token.text).join(' ');
        }

        if (next !== original) {
            const range = Range.create(
                Position.create(line.line.lineNumber, 0),
                Position.create(line.line.lineNumber, original.length)
            );
            edits.push(TextEdit.replace(range, next));
        }
    }

    return edits;
}
