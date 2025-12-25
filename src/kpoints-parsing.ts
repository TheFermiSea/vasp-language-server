import { TextDocument } from 'vscode-languageserver-textdocument';
import { Diagnostic, DiagnosticSeverity } from 'vscode-languageserver-types';

export interface KpointsData {
    comment: string;
    numKpoints: number;
    mode: string;
    grid: number[];
    shift: number[];
    isValid: boolean;
    diagnostics: Diagnostic[];
}

export function parseKpoints(document: TextDocument): KpointsData {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
    const diagnostics: Diagnostic[] = [];

    const data: KpointsData = {
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
    } else {
        data.numKpoints = numK;
    }

    // Line 3: Generation Mode
    // VASP reads the first character mainly: M, G, L, C, K
    const line3 = lines[2].trim();
    if (line3.length === 0) {
        diagnostics.push({
            severity: DiagnosticSeverity.Error,
            range: { start: { line: 2, character: 0 }, end: { line: 2, character: 0 } },
            message: 'Line 3 must specify the generation mode (e.g., Monkhorst-Pack, Gamma).'
        });
    } else {
        data.mode = line3;
    }

    // Line 4: Grid or Coordinates
    // If numKpoints == 0, Line 4 is the grid (e.g. 4 4 4)
    if (data.numKpoints === 0 && lines.length > 3) {
        const line4 = lines[3].trim();
        const tokens = line4.split(/\s+/).filter((t) => t.length > 0);
        if (tokens.length >= 3) {
            const gx = Number(tokens[0]);
            const gy = Number(tokens[1]);
            const gz = Number(tokens[2]);

            if (!Number.isInteger(gx) || !Number.isInteger(gy) || !Number.isInteger(gz)) {
                diagnostics.push({
                    severity: DiagnosticSeverity.Error,
                    range: { start: { line: 3, character: 0 }, end: { line: 3, character: line4.length } },
                    message: 'Grid values must be integers.'
                });
            } else {
                data.grid = [gx, gy, gz];
            }
        }
    }

    return data;
}
