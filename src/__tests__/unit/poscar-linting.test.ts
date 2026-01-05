import { validatePoscar } from '../../features/poscar/linting';
import { parsePoscar } from '../../features/poscar/parsing';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { DiagnosticSeverity } from 'vscode-languageserver-types';

describe('POSCAR Linter', () => {
    function validate(content: string) {
        const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
        const parsed = parsePoscar(doc);
        return validatePoscar(doc, parsed);
    }

    const validPoscar = `Example POSCAR
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
2
Direct
0.0 0.0 0.0
0.5 0.5 0.5`;

    test('passes valid POSCAR', () => {
        const diags = validate(validPoscar);
        expect(diags).toHaveLength(0);
    });

    test('validates scaling factor must be 1 or 3 values', () => {
        const invalid = `Title
1.0 2.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
2
Direct
0.0 0.0 0.0
0.5 0.5 0.5`;
        const diags = validate(invalid);
        expect(diags.length).toBeGreaterThan(0);
        expect(diags.some((d) => d.message.includes('scaling factor'))).toBeTruthy();
    });

    test('validates species names are letters only', () => {
        // Species names with numbers mixed in after valid prefix
        const invalid = `Title
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe Fe123
2 2
Direct
0.0 0.0 0.0
0.5 0.5 0.5
0.0 0.0 0.0
0.5 0.5 0.5`;
        const diags = validate(invalid);
        // "Fe123" should be flagged as invalid species name
        expect(diags.some((d) => d.message.includes('invalid'))).toBeTruthy();
    });

    test('validates atom counts are positive integers', () => {
        const invalid = `Title
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
-2
Direct
0.0 0.0 0.0`;
        const diags = validate(invalid);
        expect(diags.some((d) => d.message.includes('positive integer'))).toBeTruthy();
    });

    test('validates vector components are numbers', () => {
        const invalid = `Title
1.0
abc 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
1
Direct
0.0 0.0 0.0`;
        const diags = validate(invalid);
        expect(diags.some((d) => d.message.includes('number'))).toBeTruthy();
    });

    test('validates species/atom count mismatch', () => {
        const invalid = `Title
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe O
2
Direct
0.0 0.0 0.0
0.5 0.5 0.5`;
        const diags = validate(invalid);
        expect(diags.some((d) => d.message.includes('species'))).toBeTruthy();
    });

    test('validates selective dynamics flags', () => {
        const invalid = `Title
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
1
Selective dynamics
Direct
0.0 0.0 0.0 X Y Z`;
        const diags = validate(invalid);
        expect(diags.some((d) => d.message.includes("'T' or 'F'"))).toBeTruthy();
    });

    test('handles VASP 5 format with species names', () => {
        const vasp5 = `VASP 5 format
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe O
1 1
Direct
0.0 0.0 0.0
0.5 0.5 0.5`;
        const diags = validate(vasp5);
        expect(diags).toHaveLength(0);
    });

    test('handles per-axis scaling factors', () => {
        const threeScale = `Title
1.0 2.0 3.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
1
Direct
0.0 0.0 0.0`;
        const diags = validate(threeScale);
        expect(diags).toHaveLength(0);
    });

    test('validates negative per-axis scaling factors', () => {
        const negative = `Title
1.0 -2.0 3.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
1
Direct
0.0 0.0 0.0`;
        const diags = validate(negative);
        expect(diags.some((d) => d.message.includes('positive'))).toBeTruthy();
    });
});
