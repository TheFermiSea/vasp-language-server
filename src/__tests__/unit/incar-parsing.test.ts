import { parseIncar } from '../../incar-parsing';
import { TextDocument } from 'vscode-languageserver-textdocument';

describe('INCAR Parser', () => {
    function createDoc(content: string) {
        return TextDocument.create('file:///INCAR', 'vasp', 1, content);
    }

    test('parses simple tag=value assignments', () => {
        const doc = createDoc('ENCUT = 500\nISMEAR = 1');
        const parsed = parseIncar(doc);
        expect(parsed.statements).toHaveLength(2);
        expect(parsed.statements[0].tag.text).toBe('ENCUT');
        expect(parsed.statements[0].values[0].text).toBe('500');
    });

    test('parses multiple statements on one line', () => {
        const doc = createDoc('ENCUT = 500 ; ISMEAR = 1');
        const parsed = parseIncar(doc);
        expect(parsed.statements).toHaveLength(2);
    });

    test('ignores comments', () => {
        const doc = createDoc('ENCUT = 500 # This is a comment\n! Another comment');
        const parsed = parseIncar(doc);
        expect(parsed.statements).toHaveLength(1);
        expect(parsed.statements[0].tag.text).toBe('ENCUT');
    });

    test('handles continuations', () => {
        const doc = createDoc('SYSTEM = Long name \\\n continued');
        const parsed = parseIncar(doc);
        expect(parsed.statements).toHaveLength(1);
        expect(parsed.statements[0].values.map((t: any) => t.text).join(' ')).toContain('Long name');
        expect(parsed.statements[0].values.map((t: any) => t.text).join(' ')).toContain('continued');
    });

    test('handles empty files', () => {
        const doc = createDoc('');
        const parsed = parseIncar(doc);
        expect(parsed.statements).toHaveLength(0);
    });
});
