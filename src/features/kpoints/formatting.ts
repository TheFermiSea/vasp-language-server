import { TextDocument } from 'vscode-languageserver-textdocument';
import { Position, Range, TextEdit } from 'vscode-languageserver-types';

/**
 * Format a KPOINTS document by normalizing whitespace per line.
 *
 * @param document - LSP text document for a KPOINTS file.
 * @returns Text edits that normalize whitespace for each line.
 */
export function formatKpoints(document: TextDocument): TextEdit[] {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
    const edits: TextEdit[] = [];

    for (let i = 0; i < lines.length; i++) {
        const original = lines[i];
        if (original.length === 0) {
            continue;
        }

        let next = original;
        if (i === 0) {
            next = original.replace(/\s+$/, '');
        } else {
            next = original.trim().split(/\s+/).join(' ');
        }

        if (next !== original) {
            const range = Range.create(Position.create(i, 0), Position.create(i, original.length));
            edits.push(TextEdit.replace(range, next));
        }
    }

    return edits;
}
