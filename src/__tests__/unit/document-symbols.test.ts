import { DocumentSymbol, SymbolKind } from "vscode-languageserver-types";
import { TextDocument } from "vscode-languageserver-textdocument";
import { getIncarSymbols, getPoscarSymbols } from "../../document-symbols";

function createDoc(content: string, uri = "file:///INCAR"): TextDocument {
    return TextDocument.create(uri, "plaintext", 1, content);
}

describe('Document Symbols (Outline)', () => {
    describe('INCAR Symbols', () => {
        it('should extract tags as variables', () => {
            const content = `
                ENCUT = 500
                ISMEAR = 0
            `;
            const doc = createDoc(content);
            const symbols = getIncarSymbols(doc);

            expect(symbols).toHaveLength(2);
            expect(symbols[0].name).toBe('ENCUT');
            expect(symbols[0].kind).toBe(SymbolKind.Variable);
            expect(symbols[0].detail).toBe('500');

            expect(symbols[1].name).toBe('ISMEAR');
            expect(symbols[1].detail).toBe('0');
        });

        it('should handle ragged assignments', () => {
            const content = 'LREAL=Auto';
            const doc = createDoc(content);
            const symbols = getIncarSymbols(doc);
            expect(symbols[0].name).toBe('LREAL');
            expect(symbols[0].detail).toBe('Auto');
        });
    });

    describe('POSCAR Symbols', () => {
        it('should structure a standard POSCAR', () => {
            const content = `Fe3O4
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
Fe O
3 4
Direct
0.1 0.1 0.1
0.5 0.5 0.5
`;
            const doc = createDoc(content, "file:///POSCAR");
            const symbols = getPoscarSymbols(doc);

            expect(symbols.length).toBeGreaterThan(4);
            expect(symbols[0].name).toBe("Title/Comment");
            expect(symbols[1].name).toBe("Scaling Factor");
            expect(symbols[2].name).toBe("Lattice Vectors");

            const species = symbols.find(s => s.name === "Species Names");
            expect(species).toBeDefined();
            expect(species?.detail).toBe("Fe O");

            const coords = symbols.find(s => s.name === "Atomic Coordinates");
            expect(coords).toBeDefined();
        });
    });
});
