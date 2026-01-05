import { parsePotcar, extractBaseElement } from '../../features/potcar/parsing';
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

describe('extractBaseElement', () => {
    test('returns bare element unchanged', () => {
        expect(extractBaseElement('Fe')).toBe('Fe');
        expect(extractBaseElement('O')).toBe('O');
        expect(extractBaseElement('Ti')).toBe('Ti');
    });

    test('strips _pv suffix (p-valence)', () => {
        expect(extractBaseElement('Fe_pv')).toBe('Fe');
        expect(extractBaseElement('Ti_pv')).toBe('Ti');
    });

    test('strips _sv suffix (s+p-valence)', () => {
        expect(extractBaseElement('Ti_sv')).toBe('Ti');
        expect(extractBaseElement('W_sv')).toBe('W');
    });

    test('strips _GW suffix (GW optimized)', () => {
        expect(extractBaseElement('Fe_GW')).toBe('Fe');
        expect(extractBaseElement('O_GW')).toBe('O');
    });

    test('strips combined _sv_GW suffix', () => {
        expect(extractBaseElement('Ti_sv_GW')).toBe('Ti');
        expect(extractBaseElement('W_sv_GW')).toBe('W');
    });

    test('strips _pv_GW suffix', () => {
        expect(extractBaseElement('Fe_pv_GW')).toBe('Fe');
    });

    test('strips _d suffix (d-states)', () => {
        expect(extractBaseElement('Ge_d')).toBe('Ge');
        expect(extractBaseElement('Ga_d')).toBe('Ga');
    });

    test('strips _d_GW suffix', () => {
        expect(extractBaseElement('Ge_d_GW')).toBe('Ge');
    });

    test('strips _s suffix (soft)', () => {
        expect(extractBaseElement('O_s')).toBe('O');
        expect(extractBaseElement('N_s')).toBe('N');
    });

    test('strips _h suffix (hard)', () => {
        expect(extractBaseElement('O_h')).toBe('O');
        expect(extractBaseElement('C_h')).toBe('C');
    });

    test('strips _AE suffix (all-electron)', () => {
        expect(extractBaseElement('H_AE')).toBe('H');
        expect(extractBaseElement('He_AE')).toBe('He');
    });

    test('strips _2 and _3 suffixes (lanthanides)', () => {
        expect(extractBaseElement('Ce_2')).toBe('Ce');
        expect(extractBaseElement('Ce_3')).toBe('Ce');
        expect(extractBaseElement('Pr_3')).toBe('Pr');
    });
});

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

    // ===== Legacy Format Tests =====

    describe('Legacy Format Support', () => {
        test('Parses ultrasoft (US) format without date', () => {
            // Very old VASP format with minimal header
            const content = `   US Ti
   some other data
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Ti');
        });

        test('Parses ultrasoft (US) format with date', () => {
            const content = `   US Fe 08Apr2002
   other data
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Fe');
        });

        test('Parses concatenated ultrasoft POTCARs', () => {
            const content = `   US Ti
   End of Dataset
   US O
   End of Dataset
   US Fe 08Apr2002
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(3);
            expect(parsed.elements[0].symbol).toBe('Ti');
            expect(parsed.elements[1].symbol).toBe('O');
            expect(parsed.elements[2].symbol).toBe('Fe');
        });

        test('Parses plain PAW format (old LDA without functional)', () => {
            // PAW without _PBE/_LDA/_GGA suffix (older LDA potentials)
            const content = `   PAW Ti_sv 26Sep2005
   other data
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Ti');
        });

        test('Parses plain PAW format and strips suffix', () => {
            const content = `   PAW Fe_pv 07Sep2000
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Fe');
        });

        test('Parses TITEL tag with PAW_PBE format', () => {
            // TITEL tag format (sometimes appears without VRHFIN)
            const content = `   TITEL  = PAW_PBE Fe 06Sep2000
   LULTRA =        F
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Fe');
        });

        test('Parses TITEL tag with US format', () => {
            const content = `   TITEL = US Ti
   other data
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Ti');
        });

        test('Parses TITEL tag with PAW (plain LDA) format', () => {
            const content = `   TITEL = PAW Ti_sv 26Sep2005
   other data
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Ti');
        });

        test('Parses PAW_LDA format', () => {
            const content = `   PAW_LDA Ti 06Sep2000
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Ti');
        });

        test('Parses PAW_GGA (PW91) format', () => {
            const content = `   PAW_GGA Ti_pv 07Sep2000
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Ti');
        });

        test('Parses GW potentials with _GW suffix', () => {
            const content = `   PAW_PBE Ti_sv_GW 05Dec2013
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Ti');
        });

        test('Parses mixed legacy formats in concatenated POTCAR', () => {
            // Real-world scenario: someone might concatenate old and new potentials
            const content = `   US Ti
   End of Dataset
   PAW Fe_pv 07Sep2000
   End of Dataset
   PAW_PBE O 08Apr2002
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(3);
            expect(parsed.elements[0].symbol).toBe('Ti');
            expect(parsed.elements[1].symbol).toBe('Fe');
            expect(parsed.elements[2].symbol).toBe('O');
        });

        test('Handles elements with numeric suffixes (_2, _3)', () => {
            const content = `   PAW_PBE Ce_3 06Sep2000
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Ce');
        });

        test('Handles all-electron potentials (_AE)', () => {
            const content = `   PAW_PBE H_AE 15Jun2001
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('H');
        });

        test('VRHFIN takes priority over legacy TITEL parsing', () => {
            // When VRHFIN is present, it should be used even if other formats exist
            const content = `   PAW_PBE Fe_pv 06Sep2000
      VRHFIN = Fe: d7s1
   TITEL  = PAW_PBE Fe_pv 06Sep2000
   End of Dataset`;
            const parsed = parsePotcar(createDoc(content));
            expect(parsed.elements).toHaveLength(1);
            expect(parsed.elements[0].symbol).toBe('Fe');
            // Should use VRHFIN line number, not PAW_PBE line
            expect(parsed.elements[0].description).toContain('VRHFIN');
        });
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

    test('Falls back to CONTCAR when POSCAR missing', async () => {
        (fs.existsSync as jest.Mock).mockImplementation((testPath: string) => testPath.endsWith('CONTCAR'));
        const contcarContent = `System\n1.0\nVec1\nVec2\nVec3\nFe O\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5`;
        (fs.promises.readFile as jest.Mock).mockResolvedValue(contcarContent);

        const content = `
VRHFIN = Fe:
End
VRHFIN = O:
End
`;
        const diags = await validatePotcar(createDoc(content));
        expect(diags).toHaveLength(0);
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
