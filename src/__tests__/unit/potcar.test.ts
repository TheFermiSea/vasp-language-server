import { parsePotcar } from '../../potcar-parsing';
import { validatePotcar } from '../../potcar-linting';
import { TextDocument } from 'vscode-languageserver-textdocument';
import * as fs from 'fs';

// Mock fs for cross-file checks
jest.mock('fs');

describe('POTCAR Parser', () => {
    function createDoc(content: string) {
        return TextDocument.create('file:///POTCAR', 'vasp', 1, content);
    }

    test('Parses VRHFIN style headers', () => {
        const content = `
   PAW_PBE Fe 06Sep2000
      VRHFIN = Fe: d7s1
   TITEL  = PAW_PBE Fe 06Sep2000
   LULTRA =        F    use ultrasoft PP ?
   IUNSCR =        1    unscreen: 0-lin 1-nonlin 2-no
        `;
        const parsed = parsePotcar(createDoc(content));
        expect(parsed.elements).toHaveLength(1);
        expect(parsed.elements[0].symbol).toBe('Fe');
    });

    test('Parses concatenated POTCARs', () => {
        const content = `
   PAW_PBE Fe 06Sep2000
      VRHFIN = Fe: d7s1
   End of Potential
   PAW_PBE O 08Apr2002
      VRHFIN = O: s2p4
   End of Potential
        `;
        const parsed = parsePotcar(createDoc(content));
        expect(parsed.elements).toHaveLength(2);
        expect(parsed.elements[0].symbol).toBe('Fe');
        expect(parsed.elements[1].symbol).toBe('O');
    });

    test('Parses TITEL style (fallback)', () => {
        const content = `
   PAW_PBE Fe 06Sep2000
   DUMMY = DATA
   End of Potential
        `;
        const parsed = parsePotcar(createDoc(content));
        expect(parsed.elements).toHaveLength(1);
        expect(parsed.elements[0].symbol).toBe('Fe');
    });
});

describe('POTCAR Linter', () => {
    function createDoc(content: string) {
        return TextDocument.create('file:///test/POTCAR', 'vasp', 1, content);
    }

    test('Detects missing POSCAR', () => {
        (fs.existsSync as jest.Mock).mockReturnValue(false);
        const diags = validatePotcar(createDoc('VRHFIN = Fe:'));
        // Should have no mismatch errors if POSCAR missing
        expect(diags).toHaveLength(0);
    });

    test('Validates matching order', () => {
        (fs.existsSync as jest.Mock).mockReturnValue(true);
        // Mock POSCAR content: "Fe O\n 1.0\n..." (VASP 5 format)
        (fs.readFileSync as jest.Mock).mockReturnValue(
            `System\n1.0\nVec1\nVec2\nVec3\nFe O\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5`
        );

        const content = `
VRHFIN = Fe:
End
VRHFIN = O:
End
`;
        const diags = validatePotcar(createDoc(content));
        expect(diags).toHaveLength(0);
    });

    test('Detects mismatch', () => {
        (fs.existsSync as jest.Mock).mockReturnValue(true);
        // POSCAR says Fe O
        (fs.readFileSync as jest.Mock).mockReturnValue(`System\n1.0\n...\n...\n...\nFe O\n...`);

        // POTCAR says O Fe
        const content = `
VRHFIN = O:
End
VRHFIN = Fe:
End
`;
        const diags = validatePotcar(createDoc(content));
        expect(diags).toHaveLength(2); // Mismatch for both positions
        expect(diags[0].message).toContain("expects 'Fe'");
    });
});
