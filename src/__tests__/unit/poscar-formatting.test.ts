import { TextDocument } from 'vscode-languageserver-textdocument';
import { TextEdit } from 'vscode-languageserver-types';
import { formatPoscar } from '../../features/poscar/formatting';

function applyLineEdits(text: string, edits: TextEdit[]): string {
    const lines = text.split(/\r?\n/);
    for (const edit of edits) {
        const line = edit.range.start.line;
        lines[line] = edit.newText;
    }
    return lines.join('\n');
}

describe('POSCAR Formatting', () => {
    test('normalizes whitespace while preserving comment text', () => {
        const input = [
            'Comment line   with   spaces   ',
            '   1.0',
            ' 1  0 0',
            '0  1    0',
            '0 0   1',
            'Fe   O',
            '2  3',
            'Direct',
            '0 0 0',
            '0.5  0.5 0.5'
        ].join('\n');

        const document = TextDocument.create('file:///POSCAR', 'vasp', 1, input);
        const edits = formatPoscar(document);
        const formatted = applyLineEdits(input, edits);

        const expected = [
            'Comment line   with   spaces',
            '1.0',
            '1 0 0',
            '0 1 0',
            '0 0 1',
            'Fe O',
            '2 3',
            'Direct',
            '0 0 0',
            '0.5 0.5 0.5'
        ].join('\n');

        expect(formatted).toBe(expected);
    });
});
