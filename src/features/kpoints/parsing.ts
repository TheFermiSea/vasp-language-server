import { TextDocument } from 'vscode-languageserver-textdocument';
import { Diagnostic, DiagnosticSeverity } from 'vscode-languageserver-types';
import { splitLines } from '../../core/parser-utils';

export interface KpointsDocument {
    comment: string;
    numKpoints: number;
    mode: string;
    grid: number[];
    shift: number[];
    isValid: boolean;
    diagnostics: Diagnostic[];
}

/**
 * Parse a KPOINTS document into a structured representation with diagnostics.
 *
 * @param document - LSP text document for a KPOINTS file.
 * @returns Parsed KPOINTS document with grid/shift info and diagnostics.
 */
export function parseKpoints(document: TextDocument): KpointsDocument {
    const text = document.getText();
    const lines = splitLines(text);
    const diagnostics: Diagnostic[] = [];

    const data: KpointsDocument = {
        comment: '',
        numKpoints: 0,
        mode: '',
        grid: [],
        shift: [0, 0, 0],
        isValid: true,
        diagnostics: diagnostics
    };

    if (lines.length < 3) {
        diagnostics.push({
            severity: DiagnosticSeverity.Error,
            range: { start: { line: 0, character: 0 }, end: { line: lines.length, character: 0 } },
            message: `KPOINTS file is too short (${lines.length} line${lines.length === 1 ? '' : 's'}). A valid KPOINTS file requires at least 3 lines:\n  Line 1: Comment/description\n  Line 2: Number of k-points (0 for automatic generation)\n  Line 3: Generation mode (Gamma, Monkhorst-Pack, Auto, etc.)`
        });
        data.isValid = false;
        return data;
    }

    // Line 1: Comment
    data.comment = lines[0].trim();

    // Line 2: Number of K-Points
    const line2 = lines[1].trim();
    const firstToken = line2.split(/\s+/)[0];
    const numK = parseInt(firstToken);
    if (isNaN(numK)) {
        diagnostics.push({
            severity: DiagnosticSeverity.Error,
            range: { start: { line: 1, character: 0 }, end: { line: 1, character: line2.length } },
            message: `Line 2: Invalid k-point count '${firstToken || '(empty)'}'. Expected an integer:\n  - Use 0 for automatic mesh generation (recommended)\n  - Use a positive integer for explicit k-point specification\nExample: '0' for automatic, '4' for 4 explicit k-points`
        });
        data.isValid = false;
    } else if (numK < 0) {
        diagnostics.push({
            severity: DiagnosticSeverity.Error,
            range: { start: { line: 1, character: 0 }, end: { line: 1, character: line2.length } },
            message: `Line 2: Negative k-point count '${numK}' is invalid. Use 0 for automatic generation or a positive integer for explicit k-points.`
        });
        data.isValid = false;
    } else {
        data.numKpoints = numK;
    }

    // Line 3: Generation Mode
    const line3 = lines[2].trim().toLowerCase();
    const validModes = [
        'gamma',
        'monkhorst-pack',
        'monkhorst',
        'auto',
        'automatic',
        'line',
        'line-mode',
        'cartesian',
        'reciprocal'
    ];
    const modeFirstChar = line3.charAt(0);
    if (line3.length === 0) {
        diagnostics.push({
            severity: DiagnosticSeverity.Error,
            range: { start: { line: 2, character: 0 }, end: { line: 2, character: 0 } },
            message: `Line 3: Missing generation mode. Expected one of:\n  - 'Gamma' or 'G': Gamma-centered mesh\n  - 'Monkhorst-Pack' or 'M': Monkhorst-Pack mesh\n  - 'Auto' or 'A': Automatic length-based generation\n  - 'Line' or 'L': Line mode for band structures\n  - 'Cartesian' or 'C': Explicit k-points in Cartesian coords\n  - 'Reciprocal' or 'R': Explicit k-points in reciprocal coords`
        });
        data.isValid = false;
    } else if (!validModes.some((m) => m.startsWith(modeFirstChar))) {
        diagnostics.push({
            severity: DiagnosticSeverity.Warning,
            range: { start: { line: 2, character: 0 }, end: { line: 2, character: lines[2].trim().length } },
            message: `Line 3: Unrecognized mode '${lines[2].trim()}'. VASP will use the first character to determine mode. Common modes start with:\n  G = Gamma-centered, M = Monkhorst-Pack, A = Automatic, L = Line-mode`
        });
        data.mode = lines[2].trim();
    } else {
        data.mode = lines[2].trim();
    }

    // Line 4: Grid or Coordinates
    if (lines.length > 3) {
        const line4 = lines[3].trim();
        const tokens = line4.split(/\s+/).filter((t) => t.length > 0);

        // Mode A: Automatic length-based (Line 4 is a single float)
        if (line3.startsWith('a') && data.numKpoints === 0) {
            const length = Number(tokens[0]);
            if (tokens.length === 0 || isNaN(length)) {
                diagnostics.push({
                    severity: DiagnosticSeverity.Error,
                    range: { start: { line: 3, character: 0 }, end: { line: 3, character: line4.length } },
                    message: `Line 4: Invalid length parameter '${tokens[0] || '(empty)'}' for automatic mode. Expected a positive number representing the real-space length in Angstroms.\nExample: '20' generates a mesh with ~20 Angstrom spacing in reciprocal space.`
                });
                data.isValid = false;
            } else if (length <= 0) {
                diagnostics.push({
                    severity: DiagnosticSeverity.Error,
                    range: { start: { line: 3, character: 0 }, end: { line: 3, character: line4.length } },
                    message: `Line 4: Length must be positive (got ${length}). For automatic k-point generation, typical values are 20-50 Angstroms. Larger values = denser k-mesh.`
                });
                data.isValid = false;
            }
        }
        // Mode M/G/C: Grid (Line 4 is 3 integers)
        else if (data.numKpoints === 0) {
            if (tokens.length < 3) {
                diagnostics.push({
                    severity: DiagnosticSeverity.Error,
                    range: { start: { line: 3, character: 0 }, end: { line: 3, character: line4.length } },
                    message: `Line 4: Incomplete k-point grid. Expected 3 integers (Nx Ny Nz), found ${tokens.length} value(s): '${line4}'.\nExample: '4 4 4' for a 4x4x4 mesh, '6 6 1' for a 2D system.`
                });
                data.isValid = false;
            } else {
                const gx = Number(tokens[0]);
                const gy = Number(tokens[1]);
                const gz = Number(tokens[2]);

                const invalidValues: string[] = [];
                if (!Number.isInteger(gx)) invalidValues.push(`Nx='${tokens[0]}'`);
                if (!Number.isInteger(gy)) invalidValues.push(`Ny='${tokens[1]}'`);
                if (!Number.isInteger(gz)) invalidValues.push(`Nz='${tokens[2]}'`);

                if (invalidValues.length > 0) {
                    diagnostics.push({
                        severity: DiagnosticSeverity.Error,
                        range: { start: { line: 3, character: 0 }, end: { line: 3, character: line4.length } },
                        message: `Line 4: Grid values must be positive integers. Invalid: ${invalidValues.join(', ')}.\nExample: '4 4 4' for a uniform 4x4x4 mesh.`
                    });
                    data.isValid = false;
                } else if (gx <= 0 || gy <= 0 || gz <= 0) {
                    const zeroVals: string[] = [];
                    if (gx <= 0) zeroVals.push(`Nx=${gx}`);
                    if (gy <= 0) zeroVals.push(`Ny=${gy}`);
                    if (gz <= 0) zeroVals.push(`Nz=${gz}`);
                    diagnostics.push({
                        severity: DiagnosticSeverity.Error,
                        range: { start: { line: 3, character: 0 }, end: { line: 3, character: line4.length } },
                        message: `Line 4: Grid dimensions must be positive. Invalid: ${zeroVals.join(', ')}.\nUse at least 1 for each dimension. For 2D systems, use '1' for the non-periodic direction.`
                    });
                    data.isValid = false;
                } else {
                    data.grid = [gx, gy, gz];
                }
            }
        }
    } else if (data.numKpoints === 0) {
        diagnostics.push({
            severity: DiagnosticSeverity.Error,
            range: { start: { line: 2, character: 0 }, end: { line: 2, character: lines[2].length } },
            message: `Line 4 is required for automatic k-point generation (numKpoints=0). Expected:\n  - For Gamma/Monkhorst-Pack: 3 integers (Nx Ny Nz)\n  - For Automatic mode: 1 positive number (length)\nExample: '4 4 4' or '20'`
        });
        data.isValid = false;
    }

    // Line 5: Optional shift for automatic meshes
    if (data.numKpoints === 0 && data.grid.length === 3 && lines.length > 4) {
        const line5 = lines[4].trim();
        if (line5.length > 0) {
            const tokens = line5.split(/\s+/).filter((t) => t.length > 0);
            if (tokens.length < 3) {
                diagnostics.push({
                    severity: DiagnosticSeverity.Error,
                    range: { start: { line: 4, character: 0 }, end: { line: 4, character: line5.length } },
                    message: `Line 5: Incomplete shift vector. Expected 3 numbers (Sx Sy Sz), found ${tokens.length}: '${line5}'.\nCommon values:\n  - '0 0 0': No shift (standard for Gamma-centered)\n  - '0.5 0.5 0.5': Half-grid shift (can improve convergence)`
                });
                data.isValid = false;
            } else {
                const sx = Number(tokens[0]);
                const sy = Number(tokens[1]);
                const sz = Number(tokens[2]);
                const invalidShifts: string[] = [];
                if (Number.isNaN(sx)) invalidShifts.push(`Sx='${tokens[0]}'`);
                if (Number.isNaN(sy)) invalidShifts.push(`Sy='${tokens[1]}'`);
                if (Number.isNaN(sz)) invalidShifts.push(`Sz='${tokens[2]}'`);

                if (invalidShifts.length > 0) {
                    diagnostics.push({
                        severity: DiagnosticSeverity.Error,
                        range: { start: { line: 4, character: 0 }, end: { line: 4, character: line5.length } },
                        message: `Line 5: Invalid shift values: ${invalidShifts.join(', ')}. Shift must be numeric (typically 0 or 0.5).\nExample: '0 0 0' for no shift, '0.5 0.5 0.5' for half-grid shift.`
                    });
                    data.isValid = false;
                } else {
                    // Warn about unusual shift values
                    const unusualShifts: string[] = [];
                    if (sx < 0 || sx > 1) unusualShifts.push(`Sx=${sx}`);
                    if (sy < 0 || sy > 1) unusualShifts.push(`Sy=${sy}`);
                    if (sz < 0 || sz > 1) unusualShifts.push(`Sz=${sz}`);

                    if (unusualShifts.length > 0) {
                        diagnostics.push({
                            severity: DiagnosticSeverity.Warning,
                            range: { start: { line: 4, character: 0 }, end: { line: 4, character: line5.length } },
                            message: `Line 5: Unusual shift values: ${unusualShifts.join(', ')}. Shifts are typically between 0 and 1. Values outside this range are unusual but valid.`
                        });
                    }
                    data.shift = [sx, sy, sz];
                }
            }
        }
    }

    return data;
}
