import { TextDocument } from 'vscode-languageserver-textdocument';
import { parsePoscar } from '../../poscar-parsing';
import { parseIncar } from '../../incar-parsing';

describe('Parser Stress Tests', () => {

    test('POSCAR parsing of a massive structure (100k atoms)', () => {
        // approx 100,000 atoms * 100 bytes per line = 10MB file
        const numAtoms = 100000;
        let content = "Massive System\n 1.0\n 10.0 0.0 0.0\n 0.0 10.0 0.0\n 0.0 0.0 10.0\n C\n " + numAtoms + "\nDirect\n";

        const startTimeStr = Date.now();
        const lines: string[] = [];
        for (let i = 0; i < numAtoms; i++) {
            lines.push(" 0.12345678 0.12345678 0.12345678");
        }
        content += lines.join('\n');
        const endTimeStr = Date.now();
        console.log(`Stress: massive content generation took ${endTimeStr - startTimeStr}ms`);

        const document = TextDocument.create('file:///MASSIVE_POSCAR', 'vasp', 1, content);

        const startTime = Date.now();
        const parsed = parsePoscar(document);
        const endTime = Date.now();

        console.log(`Stress: parsing ${numAtoms} atoms took ${endTime - startTime}ms`);

        expect(parsed.lines.length).toBeGreaterThan(numAtoms);
        expect(endTime - startTime).toBeLessThan(1000); // Should parse 100k atoms in < 1s
    });

    test('INCAR parsing of a very long file (10k tags)', () => {
        const numTags = 10000;
        let content = "# Huge INCAR\n";
        for (let i = 0; i < numTags; i++) {
            content += `TAG${i} = VALUE${i} # comment ${i}\n`;
        }
        const document = TextDocument.create('file:///HUGE_INCAR', 'vasp', 1, content);

        const startTime = Date.now();
        const parsed = parseIncar(document);
        const endTime = Date.now();

        console.log(`Stress: parsing ${numTags} tags took ${endTime - startTime}ms`);

        expect(parsed.statements.length).toBe(numTags);
        expect(endTime - startTime).toBeLessThan(500); // Should parse 10k tags in < 0.5s
    });

});
