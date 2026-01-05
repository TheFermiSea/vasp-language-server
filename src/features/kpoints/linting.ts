import { KpointsDocument } from './parsing';
import { Diagnostic, DiagnosticSeverity } from 'vscode-languageserver-types';
import { createDiagnostic } from '../../utils/util';

/**
 * Validate a parsed KPOINTS document and return diagnostics.
 *
 * @param data - Parsed KPOINTS document.
 * @returns Array of diagnostics for invalid modes or grids.
 */
export function validateKpoints(data: KpointsDocument): Diagnostic[] {
    // Start with parser diagnostics
    const diagnostics: Diagnostic[] = [...data.diagnostics];

    // Validate Mode
    const mode = data.mode.toUpperCase();
    if (mode.length > 0) {
        const firstChar = mode.charAt(0);
        const validModes = ['M', 'G', 'L', 'A', 'C', 'K', 'H']; // Monkhorst, Gamma, Line, Auto, Cartesian, Reciprocal, Hexagonal
        // Note: VASP is permissive, but usually M, G, L are used for auto-grids.

        if (!validModes.includes(firstChar)) {
            diagnostics.push(
                createDiagnostic(
                    { start: { line: 2, character: 0 }, end: { line: 2, character: data.mode.length } },
                    `Unknown KPOINTS mode '${data.mode}'. Common modes: Monkhorst-Pack, Gamma, Line-mode.`,
                    DiagnosticSeverity.Warning
                )
            );
        }
    }

    // Validate Grid if auto-generation
    if (data.numKpoints === 0) {
        // Grid check
        if (data.grid.length === 3) {
            if (data.grid.some((v: number) => v <= 0)) {
                diagnostics.push(
                    createDiagnostic(
                        { start: { line: 3, character: 0 }, end: { line: 3, character: 10 } },
                        'Grid values must be positive integers.',
                        DiagnosticSeverity.Error
                    )
                );
            }
        } else if (data.mode.toUpperCase().startsWith('L')) {
            // Line-mode
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
