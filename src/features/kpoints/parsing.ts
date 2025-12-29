import { TextDocument } from 'vscode-languageserver-textdocument';
import { Diagnostic, DiagnosticSeverity } from 'vscode-languageserver-types';

export interface KpointsDocument {
    comment: string;
    numKpoints: number;
    mode: string;
    grid: number[];
    shift: number[];
    isValid: boolean;
    diagnostics: Diagnostic[];
}

export function parseKpoints(document: TextDocument): KpointsDocument {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
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
            message: 'KPOINTS file is too short. Expected at least 3 lines.'
        });
        data.isValid = false;
        return data;
    }

    // Line 1: Comment
    data.comment = lines[0].trim();

    // Line 2: Number of K-Points
    const line2 = lines[1].trim();
    const numK = parseInt(line2.split(/\s+/)[0]);
    if (isNaN(numK)) {
        diagnostics.push({
            severity: DiagnosticSeverity.Error,
            range: { start: { line: 1, character: 0 }, end: { line: 1, character: line2.length } },
            message: 'Line 2 must be an integer (0 for auto-generation, >0 for explicit).'
        });
        data.isValid = false;
    } else {
        data.numKpoints = numK;
    }

    // Line 3: Generation Mode
    const line3 = lines[2].trim().toLowerCase();
    if (line3.length === 0) {
        diagnostics.push({
            severity: DiagnosticSeverity.Error,
            range: { start: { line: 2, character: 0 }, end: { line: 2, character: 0 } },
            message: 'Line 3 must specify the generation mode (e.g., Monkhorst-Pack, Gamma, Auto).'
        });
        data.isValid = false;
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
            if (tokens.length === 0 || isNaN(length) || length <= 0) {
                diagnostics.push({
                    severity: DiagnosticSeverity.Error,
                    range: { start: { line: 3, character: 0 }, end: { line: 3, character: line4.length } },
                    message: 'Line 4 must be a positive number (length) for automatic k-point generation.'
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
                    message: 'Line 4 must define a 3-integer grid (Nx Ny Nz) for automatic meshes.'
                });
                data.isValid = false;
            } else {
                const gx = Number(tokens[0]);
                const gy = Number(tokens[1]);
                const gz = Number(tokens[2]);

                if (!Number.isInteger(gx) || !Number.isInteger(gy) || !Number.isInteger(gz)) {
                    diagnostics.push({
                        severity: DiagnosticSeverity.Error,
                        range: { start: { line: 3, character: 0 }, end: { line: 3, character: line4.length } },
                        message: 'Grid values must be integers.'
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
            message: 'Missing line 4 with grid/length definition for automatic KPOINTS.'
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
                    message: 'Shift line must contain three numbers (e.g., 0 0 0).'
                });
                data.isValid = false;
            } else {
                const sx = Number(tokens[0]);
                const sy = Number(tokens[1]);
                const sz = Number(tokens[2]);
                if ([sx, sy, sz].some((v) => Number.isNaN(v))) {
                    diagnostics.push({
                        severity: DiagnosticSeverity.Error,
                        range: { start: { line: 4, character: 0 }, end: { line: 4, character: line5.length } },
                        message: 'Shift values must be numeric.'
                    });
                    data.isValid = false;
                } else {
                    data.shift = [sx, sy, sz];
                }
            }
        }
    }

    return data;
}
