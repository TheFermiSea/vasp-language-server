import { KpointsData } from './kpoints-parsing';
import { Diagnostic, DiagnosticSeverity } from 'vscode-languageserver-types';

export function validateKpoints(data: KpointsData): Diagnostic[] {
    // Start with parser diagnostics
    const diagnostics: Diagnostic[] = [...data.diagnostics];

    // Validate Mode
    const mode = data.mode.toUpperCase();
    if (mode.length > 0) {
        const firstChar = mode.charAt(0);
        const validModes = ['M', 'G', 'L', 'A', 'C', 'K', 'H']; // Monkhorst, Gamma, Line, Auto, Cartesian, Reciprocal, Hexagonal
        // Note: VASP is permissive, but usually M, G, L are used for auto-grids.

        if (!validModes.includes(firstChar)) {
            diagnostics.push({
                severity: DiagnosticSeverity.Warning,
                range: { start: { line: 2, character: 0 }, end: { line: 2, character: data.mode.length } },
                message: `Unknown KPOINTS mode '${data.mode}'. Common modes: Monkhorst-Pack, Gamma, Line-mode.`
            });
        }
    }

    // Validate Grid if auto-generation
    if (data.numKpoints === 0) {
        // Grid check
        if (data.grid.length === 3) {
            if (data.grid.some(v => v <= 0)) {
                diagnostics.push({
                    severity: DiagnosticSeverity.Error,
                    range: { start: { line: 3, character: 0 }, end: { line: 3, character: 10 } }, // Approx range
                    message: "Grid values must be positive integers."
                });
            }
        } else if (data.mode.toUpperCase().startsWith('L')) { // Line-mode
            // Line mode doesn't use a 3-val grid on line 4 usually?
            // Actually line mode with 0 kpoints implies... wait.
            // Line-mode typically specifies intersections.
            // Usually Line-mode involves numKpoints > 0 (points per segment).
            // If numKpoints == 0, it means automatic mesh.
            // M/G require mesh.
        }
    }

    return diagnostics;
}
