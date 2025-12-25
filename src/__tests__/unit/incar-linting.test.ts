import { validateIncar } from '../../incar-linting';
import { parseIncar } from '../../incar-parsing';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { DiagnosticSeverity } from 'vscode-languageserver-types';

describe('INCAR Linter', () => {
    function validate(content: string) {
        const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
        const parsed = parseIncar(doc);
        return validateIncar(parsed);
    }

    test('validates types correctly', () => {
        const diags = validate('ENCUT = High');
        expect(diags).toHaveLength(1);
        expect(diags[0].message).toContain('Expected number');
        expect(diags[0].severity).toBe(DiagnosticSeverity.Error);
    });

    test('passes valid input', () => {
        const diags = validate('ENCUT = 500');
        expect(diags).toHaveLength(0);
    });

    test('validates boolean values', () => {
        const diags = validate('LREAL = .True.');
        expect(diags).toHaveLength(0);

        const diagsFail = validate('LREAL = Maybe');
        expect(diagsFail).toHaveLength(1);
    });

    test('detects unknown tags', () => {
        const diags = validate('UNKNOWN_TAG = 123');
        expect(diags).toHaveLength(1);
        expect(diags[0].message).toContain('Unknown tag');
    });

    test('Easter Egg: POTIM Tempo Warning', () => {
        const diags = validate('POTIM = 4.0');
        expect(diags.some((d) => d.message.includes('Presto'))).toBeTruthy();
    });
});
