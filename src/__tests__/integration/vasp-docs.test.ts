import { TextDocument } from 'vscode-languageserver-textdocument';
import { validateIncar } from '../../features/incar/linting';
import { parseIncar } from '../../features/incar/parsing';

describe('VASP Documentation Cross-Validation (Expert Mode)', () => {
    test('validates Hybrid Functional (HSE) tags (VASP 6.4.3+)', () => {
        const text = `
SYSTEM = HSE functional
LHFCALC = .TRUE.
ALGO = All
HFSCREEN = 0.2
XC = HSE06
`;
        const doc = TextDocument.create('file:///INCAR', 'vasp', 1, text);
        const parsed = parseIncar(doc);
        const diagnostics = validateIncar(parsed);

        // Should have zero diagnostics for valid HSE setup
        expect(diagnostics.filter((d) => d.severity === 1)).toHaveLength(0);
    });

    test('validates LDA+U (Hubbard) array syntax', () => {
        const text = `
SYSTEM = NiO with LDA+U
LDAU = .TRUE.
LDAUTYPE = 2
LDAUL = 2 -1
LDAUU = 5.0 0.0
LDAUJ = 0.0 0.0
`;
        const doc = TextDocument.create('file:///INCAR', 'vasp', 1, text);
        const parsed = parseIncar(doc);
        const diagnostics = validateIncar(parsed);

        expect(diagnostics.filter((d) => d.severity === 1)).toHaveLength(0);
    });

    test('handles VASP repetition notation (N*val)', () => {
        const text = `
MAGMOM = 10*0.0 5*1.2
ENCUT = 10*500
`;
        const doc = TextDocument.create('file:///INCAR', 'vasp', 1, text);
        const parsed = parseIncar(doc);
        const diagnostics = validateIncar(parsed);

        // 10*0.0 is valid for array, but 10*500 is technically valid as VASP parses first item
        expect(diagnostics.filter((d) => d.severity === 1)).toHaveLength(0);
    });

    test('validates van der Waals (IVDW) settings', () => {
        const text = `
IVDW = 11
LVDW_EWALD = .TRUE.
`;
        const doc = TextDocument.create('file:///INCAR', 'vasp', 1, text);
        const parsed = parseIncar(doc);
        const diagnostics = validateIncar(parsed);

        // LVDW_EWALD is a valid sub-tag (though not in my common list yet, it should warn if unknown)
        const unknownTags = diagnostics.filter((d) => d.message.includes('Unknown tag'));
        expect(unknownTags.some((d) => d.message.includes('LVDW_EWALD'))).toBe(true);

        // IVDW = 11 should be valid
        const ivdwErrors = diagnostics.filter((d) => d.message.includes('IVDW'));
        expect(ivdwErrors).toHaveLength(0);
    });

    test('handles negative values in arrays (e.g. EDIFFG or POTIM)', () => {
        const text = `
EDIFFG = -0.01
POTIM = 0.5
`;
        const doc = TextDocument.create('file:///INCAR', 'vasp', 1, text);
        const parsed = parseIncar(doc);
        const diagnostics = validateIncar(parsed);
        expect(diagnostics.filter((d) => d.severity === 1)).toHaveLength(0);
    });

    test('handles leading and trailing dots for floats', () => {
        const text = `
ENCUT = 500.
SIGMA = .05
`;
        const doc = TextDocument.create('file:///INCAR', 'vasp', 1, text);
        const parsed = parseIncar(doc);
        const diagnostics = validateIncar(parsed);
        expect(diagnostics.filter((d) => d.severity === 1)).toHaveLength(0);
    });

    test('handles Fortran D-notation for floats', () => {
        const text = `
ENCUT = 5.0D+2
SIGMA = -1.0d-2
`;
        const doc = TextDocument.create('file:///INCAR', 'vasp', 1, text);
        const parsed = parseIncar(doc);
        const diagnostics = validateIncar(parsed);
        expect(diagnostics.filter((d) => d.severity === 1)).toHaveLength(0);
    });
});

import { parseKpoints } from '../../features/kpoints/parsing';

describe('VASP Documentation KPOINTS Validation', () => {
    test('validates Automatic length-based KPOINTS (Mode A)', () => {
        const text = `Automatic mesh
0
Auto
40
`;
        const doc = TextDocument.create('file:///KPOINTS', 'vasp', 1, text);
        const parsed = parseKpoints(doc);
        expect(parsed.diagnostics.filter((d) => d.severity === 1)).toHaveLength(0);
    });

    test('validates Monkhorst-Pack grid (Mode M)', () => {
        const text = `Regular mesh
0
Monkhorst-Pack
4 4 4
0 0 0
`;
        const doc = TextDocument.create('file:///KPOINTS', 'vasp', 1, text);
        const parsed = parseKpoints(doc);
        expect(parsed.diagnostics.filter((d) => d.severity === 1)).toHaveLength(0);
        expect(parsed.grid).toEqual([4, 4, 4]);
    });

    test('fails on invalid grid integers', () => {
        const text = `Broken mesh
0
Gamma
4.5 4 4
`;
        const doc = TextDocument.create('file:///KPOINTS', 'vasp', 1, text);
        const parsed = parseKpoints(doc);
        // Updated error message now says "positive integers"
        expect(parsed.diagnostics.some((d) => d.message.match(/must be.*integers/i))).toBe(true);
    });

    test('fails on invalid auto length', () => {
        const text = `Broken auto
0
Automatic
-5.0
`;
        const doc = TextDocument.create('file:///KPOINTS', 'vasp', 1, text);
        const parsed = parseKpoints(doc);
        // Updated error message now says "must be positive"
        expect(parsed.diagnostics.some((d) => d.message.match(/must be positive|positive/i))).toBe(true);
    });
});
