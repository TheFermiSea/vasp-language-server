import { getFoldingRanges } from '../../features/folding';
import { TextDocument } from 'vscode-languageserver-textdocument';

describe('Folding Ranges', () => {
    it('should fold Lattice Vectors and Coordinates in POSCAR', () => {
        // Mock a standard POSCAR
        const content = [
            'Si', // 0: Title
            '1.0', // 1: Scaling
            '5.4 0.0 0.0', // 2: Lattice Start
            '0.0 5.4 0.0', // 3
            '0.0 0.0 5.4', // 4: Lattice End
            'Si', // 5: Species
            '2', // 6: Counts
            'D', // 7: Coord Sys
            '0.0 0.0 0.0', // 8: Atom 1
            '0.25 0.25 0.25' // 9: Atom 2
        ].join('\n');

        const document = TextDocument.create('file:///test/POSCAR', 'vasp', 1, content);
        const ranges = getFoldingRanges(document);

        // Expect 2 ranges: Lattice (2-4) and Coords (8-9)
        expect(ranges.length).toBeGreaterThanOrEqual(2);

        // Lattice
        const lattice = ranges.find((r) => r.startLine === 2 && r.endLine === 4);
        expect(lattice).toBeDefined();

        // Coordinates
        const coords = ranges.find((r) => r.startLine === 8 && r.endLine === 9);
        expect(coords).toBeDefined();
    });

    it('should handle small POSCARs gracefully', () => {
        const content = 'Title\n1.0\n';
        const document = TextDocument.create('file:///test/POSCAR', 'vasp', 1, content);
        const ranges = getFoldingRanges(document);
        expect(ranges).toHaveLength(0);
    });

    it('should return empty for unknown file types', () => {
        const document = TextDocument.create('file:///test/UNKNOWN', 'vasp', 1, 'content');
        const ranges = getFoldingRanges(document);
        expect(ranges).toHaveLength(0);
    });
});
