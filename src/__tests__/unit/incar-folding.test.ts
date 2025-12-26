import { TextDocument } from 'vscode-languageserver-textdocument';
import { getIncarFoldingRanges } from '../../features/incar/folding';
import { FoldingRangeKind } from 'vscode-languageserver-types';

describe('INCAR Folding', () => {
    it('should fold sections based on headers', () => {
        const text = `# --- Electronic ---
ENCUT = 500
EDIFF = 1E-5

# --- Ionic ---
NSW = 100
ISIF = 3`;
        const doc = TextDocument.create('file:///INCAR', 'vasp', 1, text);
        const ranges = getIncarFoldingRanges(doc);

        expect(ranges).toHaveLength(2);
        expect(ranges[0].startLine).toBe(0);
        expect(ranges[0].endLine).toBe(3); // Ends before next header
        expect(ranges[0].kind).toBe(FoldingRangeKind.Region);

        expect(ranges[1].startLine).toBe(4);
        expect(ranges[1].endLine).toBe(6); // Ends at EOF
    });

    it('should handle nested or complex files gracefully', () => {
        const text = `SYSTEM = Test
        
# Section 1
TAG = 1`;
        const doc = TextDocument.create('file:///INCAR', 'vasp', 1, text);
        const ranges = getIncarFoldingRanges(doc);

        expect(ranges).toHaveLength(1);
        expect(ranges[0].startLine).toBe(2);
        expect(ranges[0].endLine).toBe(3);
    });
});
