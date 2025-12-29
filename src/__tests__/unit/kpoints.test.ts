import { parseKpoints } from '../../features/kpoints/parsing';
import { validateKpoints } from '../../features/kpoints/linting';
import { TextDocument } from 'vscode-languageserver-textdocument';

describe('KPOINTS Parser', () => {
    function createDoc(content: string) {
        return TextDocument.create('file:///KPOINTS', 'vasp', 1, content);
    }

    test('Parses valid Monkhorst-Pack', () => {
        const content = `Automatic Generation
0
Monkhorst-Pack
4 4 4
0 0 0`;
        const parsed = parseKpoints(createDoc(content));
        expect(parsed.isValid).toBe(true);
        expect(parsed.numKpoints).toBe(0);
        expect(parsed.mode).toContain('Monkhorst');
        expect(parsed.grid).toEqual([4, 4, 4]);
    });

    test('Parses shift line for automatic grids', () => {
        const content = `Automatic Generation
0
Gamma
4 4 4
0.5 0 0`;
        const parsed = parseKpoints(createDoc(content));
        expect(parsed.isValid).toBe(true);
        expect(parsed.shift).toEqual([0.5, 0, 0]);
    });

    test('Fails on short file', () => {
        const content = `Header
0`;
        const parsed = parseKpoints(createDoc(content));
        expect(parsed.isValid).toBe(false);
    });

    test('Parses Explicit Kpoints (Line 3 is not mode)', () => {
        const content = `Explicit
10
Cartesian
0.1 0.2 0.3 1.0
...`;
        const parsed = parseKpoints(createDoc(content));
        expect(parsed.numKpoints).toBe(10);
        expect(parsed.mode).toBe('Cartesian');
    });
});

describe('KPOINTS Linter', () => {
    function createDoc(content: string) {
        return TextDocument.create('file:///KPOINTS', 'vasp', 1, content);
    }

    test('Validates Mode String', () => {
        const content = `Header
0
ZebraMode
4 4 4`;
        const parsed = parseKpoints(createDoc(content));
        const diags = validateKpoints(parsed);
        expect(diags).toHaveLength(1);
        expect(diags[0].message).toContain('Unknown KPOINTS mode');
    });

    test('Validates Integer Grid', () => {
        const content = `Header
0
Gamma
4.5 4 4`;
        const parsed = parseKpoints(createDoc(content));
        const diags = validateKpoints(parsed);
        // Parser might flag this as NaN or linter checks it
        // The parser expects integers for the grid logic I wrote
        expect(diags.length).toBeGreaterThan(0);
        expect(diags[0].message).toContain('must be integers');
    });
});
