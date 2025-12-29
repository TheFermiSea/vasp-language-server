import { parsePotcar } from '../../features/potcar/parsing';
import { validatePotcar } from '../../features/potcar/linting';
import { TextDocument } from 'vscode-languageserver-textdocument';
import * as fs from 'fs';

// Mock fs for cross-file checks
jest.mock('fs', () => ({
    existsSync: jest.fn(),
    promises: {
        readFile: jest.fn()
    },
    readFileSync: jest.fn() // For some sync path if any left
}));

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
});

describe('POTCAR Linter', () => {
    function createDoc(content: string) {
        return TextDocument.create('file:///test/POTCAR', 'vasp', 1, content);
    }

    test('Detects missing POSCAR', async () => {
        (fs.existsSync as jest.Mock).mockReturnValue(false);
        const diags = await validatePotcar(createDoc('VRHFIN = Fe:'));
        expect(diags).toHaveLength(0);
    });

    test('Validates matching order', async () => {
        (fs.existsSync as jest.Mock).mockReturnValue(true);
        const poscarContent = `System\n1.0\nVec1\nVec2\nVec3\nFe O\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5`;
        (fs.promises.readFile as jest.Mock).mockResolvedValue(poscarContent);

        const content = `
VRHFIN = Fe:
End
VRHFIN = O:
End
`;
        const diags = await validatePotcar(createDoc(content));
        expect(diags).toHaveLength(0);
    });

    test('Detects mismatch', async () => {
        (fs.existsSync as jest.Mock).mockReturnValue(true);
        // POSCAR says Fe O
        (fs.promises.readFile as jest.Mock).mockResolvedValue(`System\n1.0\n...\n...\n...\nFe O\n...`);

        // POTCAR says O Fe
        const content = `
VRHFIN = O:
End
VRHFIN = Fe:
End
`;
        const diags = await validatePotcar(createDoc(content));
        expect(diags).toHaveLength(2);
        expect(diags[0].message).toContain("expects 'Fe'");
    });

    test('Detects no elements in junk file', async () => {
        const content = 'total garbage';
        const doc = createDoc(content);
        const diags = await validatePotcar(doc);
        expect(diags.length).toBe(1);
        expect(diags[0].message).toContain('No elements');
    });

    test('Detects missing element in POTCAR', async () => {
        const content = 'POTCAR content\n VRHFIN = Fe:\n VRHFIN = O:';
        const doc = createDoc(content);

        (fs.promises.readFile as jest.Mock).mockResolvedValue('POSCAR content\n1.0\nVec1\nVec2\nVec3\nFe O H\n 1 1 1');
        (fs.existsSync as jest.Mock).mockReturnValue(true);

        const diags = await validatePotcar(doc);
        expect(diags.length).toBe(1);
        expect(diags[0].message).toContain('Missing potential');
    });
});
