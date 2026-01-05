import { validatePoscar } from '../../features/poscar/linting';
import { parsePoscar } from '../../features/poscar/parsing';
import { validateIncar } from '../../features/incar/linting';
import { parseIncar } from '../../features/incar/parsing';
import { TextDocument } from 'vscode-languageserver-textdocument';

describe('Error Scenarios - Malformed Input', () => {
    describe('Empty and Minimal Files', () => {
        test('handles empty POSCAR', () => {
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, '');
            const parsed = parsePoscar(doc);
            const diags = validatePoscar(doc, parsed);
            // Should not crash, may have diagnostics
            expect(parsed).toBeDefined();
            expect(Array.isArray(diags)).toBe(true);
        });

        test('handles empty INCAR', () => {
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, '');
            const parsed = parseIncar(doc);
            const diags = validateIncar(parsed);
            expect(parsed).toBeDefined();
            expect(diags).toHaveLength(0);
        });

        test('handles single line POSCAR', () => {
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, 'Title only');
            const parsed = parsePoscar(doc);
            const diags = validatePoscar(doc, parsed);
            expect(parsed).toBeDefined();
            expect(parsed.lines.length).toBe(1);
        });

        test('handles whitespace-only INCAR', () => {
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, '   \n\t\n   ');
            const parsed = parseIncar(doc);
            const diags = validateIncar(parsed);
            expect(diags).toHaveLength(0);
        });
    });

    describe('Truncated Files', () => {
        test('handles POSCAR with only title and scaling', () => {
            const content = `Title
1.0`;
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
            const parsed = parsePoscar(doc);
            expect(parsed).toBeDefined();
        });

        test('handles POSCAR with incomplete lattice vectors', () => {
            const content = `Title
1.0
3.0 0.0 0.0`;
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
            const parsed = parsePoscar(doc);
            expect(parsed).toBeDefined();
        });

        test('handles POSCAR without positions', () => {
            const content = `Title
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
2
Direct`;
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
            const parsed = parsePoscar(doc);
            const diags = validatePoscar(doc, parsed);
            expect(parsed).toBeDefined();
        });
    });

    describe('Malformed Values', () => {
        test('handles INCAR with extremely long lines', () => {
            const longValue = 'A'.repeat(10000);
            const content = `SYSTEM = ${longValue}`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            const diags = validateIncar(parsed);
            expect(parsed).toBeDefined();
        });

        test('handles POSCAR with special characters', () => {
            const content = `Title with @#$%^&*()
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
1
Direct
0.0 0.0 0.0`;
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
            const parsed = parsePoscar(doc);
            expect(parsed).toBeDefined();
        });

        test('handles INCAR with Unicode characters', () => {
            const content = `SYSTEM = Test with émojis 🧪`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            expect(parsed).toBeDefined();
        });

        test('handles scientific notation edge cases', () => {
            const content = `Title
1.0E-300
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
1
Direct
1.0D+308 -1.0E-308 0.0`;
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
            const parsed = parsePoscar(doc);
            expect(parsed).toBeDefined();
        });
    });

    describe('Unusual Formatting', () => {
        test('handles INCAR with mixed line endings', () => {
            const content = `ENCUT = 500\r\nISMEAR = 0\rSIGMA = 0.1\nEDIFF = 1E-6`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            expect(parsed.statements.length).toBeGreaterThan(0);
        });

        test('handles INCAR with tabs and spaces mixed', () => {
            const content = `ENCUT\t=\t500\n  ISMEAR  =  0`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            const diags = validateIncar(parsed);
            expect(parsed).toBeDefined();
            expect(diags).toHaveLength(0);
        });

        test('handles POSCAR with trailing whitespace', () => {
            const content = `Title
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
1
Direct
0.0 0.0 0.0   `;
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
            const parsed = parsePoscar(doc);
            const diags = validatePoscar(doc, parsed);
            expect(parsed).toBeDefined();
        });

        test('handles INCAR with line continuation', () => {
            const content = `MAGMOM = 1 2 3 \\
4 5 6`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            expect(parsed).toBeDefined();
        });
    });
});
