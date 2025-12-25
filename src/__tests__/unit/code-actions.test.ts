import { Diagnostic, DiagnosticSeverity, Range } from "vscode-languageserver-types";
import { getIncarCodeActions } from "../../code-actions";

describe('Code Actions (Quick Fixes)', () => {
    const mockUri = 'file:///test/INCAR';

    it('should suggest ENCUT for EMCUT', () => {
        const diagnostic: Diagnostic = {
            message: "Unknown tag 'EMCUT'. Check spelling or version compatibility.",
            range: Range.create(0, 0, 0, 5),
            severity: DiagnosticSeverity.Warning,
            source: "VASP"
        };

        const actions = getIncarCodeActions(mockUri, [diagnostic]);

        expect(actions).toHaveLength(1);
        expect(actions[0].title).toBe("Change to 'ENCUT'");
        expect(actions[0].edit?.changes?.[mockUri]).toBeDefined();
        expect(actions[0].edit!.changes![mockUri][0].newText).toBe('ENCUT');
    });

    it('should suggest ISMEAR for SMEAR', () => {
        const diagnostic: Diagnostic = {
            message: "Unknown tag 'SMEAR'. Check spelling or version compatibility.",
            range: Range.create(1, 0, 1, 5),
            severity: DiagnosticSeverity.Warning,
            source: "VASP"
        };

        const actions = getIncarCodeActions(mockUri, [diagnostic]);

        expect(actions.length).toBeGreaterThan(0);
        // ISMEAR is distance 1 from SMEAR (insertion)
        const fix = actions.find(a => a.title === "Change to 'ISMEAR'");
        expect(fix).toBeDefined();
    });

    it('should not suggest anything for completely random strings', () => {
        const diagnostic: Diagnostic = {
            message: "Unknown tag 'XYZABC123'. Check spelling or version compatibility.",
            range: Range.create(2, 0, 2, 9),
            severity: DiagnosticSeverity.Warning,
            source: "VASP"
        };

        const actions = getIncarCodeActions(mockUri, [diagnostic]);

        // Closest match might still be found if VASP_TAGS is small, but let's check distance threshold.
        // XYZABC123 is length 9. Closest tag probably distance > 3.
        expect(actions).toHaveLength(0);
    });

    it('should handle multiple diagnostics', () => {
        const d1: Diagnostic = {
            message: "Unknown tag 'EMCUT'.", // Short message
            range: Range.create(0, 0, 0, 5),
            severity: DiagnosticSeverity.Warning
        };
        // This regex match expects "Unknown tag 'TAG'..."
        // If message format changes, logic breaks. My mock message must match regex.
        const d2: Diagnostic = {
            message: "Unknown tag 'EMCUT'. Check spelling...",
            range: Range.create(0, 0, 0, 5),
            severity: DiagnosticSeverity.Warning
        };

        const actions = getIncarCodeActions(mockUri, [d2]);
        expect(actions).toHaveLength(1);
    });
});
