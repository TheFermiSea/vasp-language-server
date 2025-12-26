import { formatIncar } from '../../features/incar/formatting';
import { TextDocument } from 'vscode-languageserver-textdocument';

describe('INCAR Formatter', () => {
    function createDoc(content: string) {
        return TextDocument.create('file:///INCAR', 'vasp', 1, content);
    }

    test('Aligns simple assignments', () => {
        const content = `ENCUT = 500
ISMEAR = 0
LREAL = Auto`;
        const expected = `ENCUT  = 500
ISMEAR = 0
LREAL  = Auto`;

        const doc = createDoc(content);
        const edits = formatIncar(doc);
        const newText = TextDocument.applyEdits(doc, edits);

        expect(newText).toBe(expected);
    });

    test('Aligns mixed with comments', () => {
        const content = `ENCUT=500
# Comment here
  SIGMA =  0.05 ! Broadening`;
        // Max length is ENCUT (5) and SIGMA (5). No, SIGMA has leading spaces?
        // My parser uses line.substring(0, eqIndex).trim().
        // So keys are "ENCUT", "SIGMA". Max len 5.
        // Padded key: "ENCUT", "SIGMA" (both length 5).
        // Wait, "ENCUT" len 5. "SIGMA" len 5.
        // Output format: KEY   = VALUE
        // "ENCUT = 500"
        // "SIGMA = 0.05 ! Broadening" (Note: trims value too?)
        // The implementation does: const value = line.substring(eqIndex + 1).trim();

        const expected = `ENCUT = 500
# Comment here
SIGMA = 0.05 ! Broadening`;

        const doc = createDoc(content);
        const edits = formatIncar(doc);
        const newText = TextDocument.applyEdits(doc, edits);

        expect(newText).toBe(expected);
    });

    test('Handles Ragged', () => {
        const content = `A=1
LONG_TAG=2`;
        // Max key: LONG_TAG (8)
        // A needs padding: A (1) + 7 spaces = 8.
        const expected = `A        = 1
LONG_TAG = 2`;
        const doc = createDoc(content);
        const edits = formatIncar(doc);
        const newText = TextDocument.applyEdits(doc, edits);

        expect(newText).toBe(expected);
    });
});
