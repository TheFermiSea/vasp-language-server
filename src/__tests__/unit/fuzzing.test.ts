import { TextDocument } from 'vscode-languageserver-textdocument';
import { parsePoscar } from '../../poscar-parsing';
import { parseIncar } from '../../incar-parsing';
import { parseKpoints } from '../../kpoints-parsing';

describe('Fuzzing-lite (Robustness Tests)', () => {

    function generateRandomGarbage(length: number): string {
        const chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 \n\t!@#$%^&*()_+-=[]{}|;:,.<>?/\\';
        let result = '';
        for (let i = 0; i < length; i++) {
            result += chars.charAt(Math.floor(Math.random() * chars.length));
        }
        return result;
    }

    test('Parsers should not crash on random garbage', () => {
        for (let i = 0; i < 100; i++) {
            const garbage = generateRandomGarbage(1000);
            const doc = TextDocument.create(`file:///fuzz_${i}`, 'vasp', 1, garbage);

            // Just ensure they don't throw
            expect(() => parseIncar(doc)).not.toThrow();
            expect(() => parsePoscar(doc)).not.toThrow();
            expect(() => parseKpoints(doc)).not.toThrow();
        }
    });

    test('Parsers should handle extremely long single lines', () => {
        const longLine = 'A'.repeat(1000000); // 1MB line
        const doc = TextDocument.create('file:///LONG_LINE', 'vasp', 1, longLine);

        expect(() => parseIncar(doc)).not.toThrow();
        expect(() => parsePoscar(doc)).not.toThrow();
        expect(() => parseKpoints(doc)).not.toThrow();
    });

    test('Parsers should handle deeply nested? (N/A for flat VASP but good to test recursively if any)', () => {
        // VASP is flat, but let's test many empty lines
        const content = '\n'.repeat(10000);
        const doc = TextDocument.create('file:///EMPTY_LINES', 'vasp', 1, content);

        expect(() => parseIncar(doc)).not.toThrow();
        expect(() => parsePoscar(doc)).not.toThrow();
        expect(() => parseKpoints(doc)).not.toThrow();
    });

});
