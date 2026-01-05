import { TextDocument } from 'vscode-languageserver-textdocument';
import { parseIncar } from '../../features/incar/parsing';
import { parsePoscar } from '../../features/poscar/parsing';
import { parseKpoints } from '../../features/kpoints/parsing';

describe('Parser Fuzzing - Edge Inputs', () => {
    test('parsers handle binary-like content', () => {
        const bytes = Array.from({ length: 256 }, (_, i) => String.fromCharCode(i)).join('');
        const doc = TextDocument.create('file:///BINARY', 'vasp', 1, bytes);

        expect(() => parseIncar(doc)).not.toThrow();
        expect(() => parsePoscar(doc)).not.toThrow();
        expect(() => parseKpoints(doc)).not.toThrow();
    });

    test('parsers handle BOM and zero-width characters', () => {
        const content = `\uFEFFTITLE\u200B\n0\nGamma\n4 4 4\n0 0 0`;
        const doc = TextDocument.create('file:///UNICODE', 'vasp', 1, content);

        expect(() => parseIncar(doc)).not.toThrow();
        expect(() => parsePoscar(doc)).not.toThrow();
        expect(() => parseKpoints(doc)).not.toThrow();
    });
});
