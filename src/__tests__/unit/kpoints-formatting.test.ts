import { TextDocument } from 'vscode-languageserver-textdocument';
import { TextEdit } from 'vscode-languageserver-types';
import { formatKpoints } from '../../features/kpoints/formatting';

function applyLineEdits(text: string, edits: TextEdit[]): string {
    const lines = text.split(/\r?\n/);
    for (const edit of edits) {
        const line = edit.range.start.line;
        lines[line] = edit.newText;
    }
    return lines.join('\n');
}

describe('KPOINTS Formatting', () => {
    test('normalizes whitespace for data lines', () => {
        const input = 'KPOINTS   comment   line   \n 0\n  Monkhorst-Pack \n  4   4  4\n  0   0  0';

        const document = TextDocument.create('file:///KPOINTS', 'vasp', 1, input);
        const edits = formatKpoints(document);
        const formatted = applyLineEdits(input, edits);

        const expected = 'KPOINTS   comment   line\n0\nMonkhorst-Pack\n4 4 4\n0 0 0';

        expect(formatted).toBe(expected);
    });
});
