import { validatePoscar } from '../../features/poscar/linting';
import { parsePoscar } from '../../features/poscar/parsing';
import { validateIncar } from '../../features/incar/linting';
import { parseIncar } from '../../features/incar/parsing';
import { validateKpoints } from '../../features/kpoints/linting';
import { parseKpoints } from '../../features/kpoints/parsing';
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

        test('handles INCAR with only comments', () => {
            const content = `# comment\n! another comment`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            const diags = validateIncar(parsed);
            expect(parsed).toBeDefined();
            expect(diags).toHaveLength(0);
        });

        test('handles POSCAR with only comments', () => {
            const content = `# comment\n! another comment`;
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
            const parsed = parsePoscar(doc);
            const diags = validatePoscar(doc, parsed);
            expect(parsed).toBeDefined();
            expect(Array.isArray(diags)).toBe(true);
        });

        test('handles single line POSCAR', () => {
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, 'Title only');
            const parsed = parsePoscar(doc);
            // Parser now returns early for files < 7 lines with diagnostics
            expect(parsed).toBeDefined();
            // Short files now produce diagnostics about being too short
            expect(parsed.diagnostics.length).toBeGreaterThan(0);
            expect(parsed.diagnostics[0].message).toContain('too short');
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
            expect(Array.isArray(diags)).toBe(true);
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
            expect(Array.isArray(diags)).toBe(true);
        });

        test('handles INCAR with unclosed multi-line string', () => {
            const content = `SYSTEM = "Unclosed string\nENCUT = 400`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            const diags = validateIncar(parsed);
            expect(parsed).toBeDefined();
            expect(diags.some((d) => d.message.includes('Unclosed string'))).toBeTruthy();
        });

        test('handles INCAR with very long tag names', () => {
            const tag = 'A'.repeat(200);
            const content = `${tag} = 1`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            const diags = validateIncar(parsed);
            expect(parsed).toBeDefined();
            expect(Array.isArray(diags)).toBe(true);
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

        test('handles POSCAR with NaN and Infinity coordinates', () => {
            const content = `Title
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
1
Direct
NaN Infinity 0.0`;
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
            const parsed = parsePoscar(doc);
            const diags = validatePoscar(doc, parsed);
            expect(parsed).toBeDefined();
            expect(diags.length).toBeGreaterThan(0);
        });

        test('handles INCAR with Unicode characters', () => {
            const content = `SYSTEM = Test with émojis 🧪`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            expect(parsed).toBeDefined();
        });

        test('handles INCAR with BOM/RTL/zero-width characters', () => {
            const content = `\uFEFFSYSTEM = value \u200B \u202E`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            const diags = validateIncar(parsed);
            expect(parsed).toBeDefined();
            expect(Array.isArray(diags)).toBe(true);
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
            expect(Array.isArray(diags)).toBe(true);
        });

        test('handles POSCAR with incomplete selective dynamics flags', () => {
            const content = `Title
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe
1
Selective dynamics
Direct
0.0 0.0 0.0 T F`;
            const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
            const parsed = parsePoscar(doc);
            const diags = validatePoscar(doc, parsed);
            expect(parsed).toBeDefined();
            expect(diags.some((d) => d.message.includes('selective-dynamics'))).toBeTruthy();
        });

        test('handles INCAR with line continuation', () => {
            const content = `MAGMOM = 1 2 3 \\
4 5 6`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            expect(parsed).toBeDefined();
        });

        test('handles INCAR with deep line continuations', () => {
            const content = `MAGMOM = 1 2 3 \\
4 5 6 \\
7 8 9`;
            const doc = TextDocument.create('file:///INCAR', 'vasp', 1, content);
            const parsed = parseIncar(doc);
            expect(parsed).toBeDefined();
        });
    });

    describe('KPOINTS edge cases', () => {
        test('handles negative k-point counts', () => {
            const content = `Comment\n-4\nGamma\n4 4 4\n0 0 0`;
            const doc = TextDocument.create('file:///KPOINTS', 'vasp', 1, content);
            const parsed = parseKpoints(doc);
            const diags = validateKpoints(parsed);
            expect(parsed).toBeDefined();
            expect(diags.length).toBeGreaterThan(0);
        });

        test('handles invalid generation modes', () => {
            const content = `Comment\n0\nZ\n4 4 4\n0 0 0`;
            const doc = TextDocument.create('file:///KPOINTS', 'vasp', 1, content);
            const parsed = parseKpoints(doc);
            const diags = validateKpoints(parsed);
            expect(parsed).toBeDefined();
            expect(diags.length).toBeGreaterThan(0);
        });
    });
});
