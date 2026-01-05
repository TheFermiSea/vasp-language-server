import { TextDocument } from 'vscode-languageserver-textdocument';
import { Diagnostic, DiagnosticSeverity } from 'vscode-languageserver-types';
import { parsePotcar } from './parsing';
import { fileURLToPath } from 'url';
import * as fs from 'fs';
import * as path from 'path';
import { logger } from '../../utils/logger';

/**
 * Validates POTCAR consistency with POSCAR.
 */
export async function validatePotcar(document: TextDocument): Promise<Diagnostic[]> {
    const diagnostics: Diagnostic[] = [];
    const parsedPot = parsePotcar(document);

    if (parsedPot.elements.length === 0) {
        diagnostics.push({
            severity: DiagnosticSeverity.Warning,
            range: { start: { line: 0, character: 0 }, end: { line: 0, character: 10 } },
            message: 'No elements detected in POTCAR. Is this a valid VASP POTCAR?'
        });
        return diagnostics;
    }

    // Attempt to find POSCAR in the same directory for cross-validation
    let filePath: string;
    try {
        filePath = fileURLToPath(document.uri);
    } catch {
        // Non-file URIs (untitled documents, remote files) cannot be cross-validated
        // Skip POSCAR validation gracefully rather than using brittle manual conversion
        return diagnostics;
    }

    const dir = path.dirname(filePath);
    const poscarPath = path.join(dir, 'POSCAR');

    if (fs.existsSync(poscarPath)) {
        try {
            const poscarContent = await fs.promises.readFile(poscarPath, 'utf-8');
            const poscarElements = getPoscarElements(poscarContent);

            if (poscarElements.length === 0) {
                // POSCAR might be implicit (old format) or parse failed
                // Skip check if we can't be sure
            } else {
                // Compare
                // 1. Count mismatch
                const maxLen = Math.max(poscarElements.length, parsedPot.elements.length);
                for (let i = 0; i < maxLen; i++) {
                    const pos = poscarElements[i];
                    const pot = parsedPot.elements[i];

                    if (!pos) {
                        // POTCAR has more elements than POSCAR
                        diagnostics.push({
                            severity: DiagnosticSeverity.Warning,
                            range: {
                                start: { line: pot.line, character: 0 },
                                end: { line: pot.line, character: pot.description.length }
                            },
                            message: `Extra potential for '${pot.symbol}' found in POTCAR but not in POSCAR.`
                        });
                    } else if (!pot) {
                        // POSCAR has more elements
                        // We can't mark strict line in POTCAR as it ends.
                        // Add error to the last POTCAR element
                        const last = parsedPot.elements[parsedPot.elements.length - 1];
                        diagnostics.push({
                            severity: DiagnosticSeverity.Error,
                            range: {
                                start: { line: last.line, character: 0 },
                                end: { line: last.line, character: last.description.length }
                            },
                            message: `Missing potential for POSCAR element '${pos}'. POTCAR ends prematurely.`
                        });
                    } else {
                        // Mismatch
                        // Handle "Fe_pv" vs "Fe": usually we just check if POSCAR symbol is a prefix or exact match
                        // Actually VASP requires strict ordering.
                        // POSCAR: "Fe O", POTCAR must be Fe then O.
                        // But POTCAR might be "Fe_pv".
                        if (!pot.symbol.startsWith(pos) && !pos.startsWith(pot.symbol)) {
                            diagnostics.push({
                                severity: DiagnosticSeverity.Error,
                                range: {
                                    start: { line: pot.line, character: 0 },
                                    end: { line: pot.line, character: pot.description.length }
                                },
                                message: `Mismatch: POSCAR expects '${pos}', but found '${pot.symbol}' in POTCAR.`
                            });
                        }
                    }
                }
            }
        } catch (error) {
            logger.warn(
                `Could not read POSCAR for cross-validation: ${error instanceof Error ? error.message : 'unknown error'}`
            );
        }
    }

    return diagnostics;
}

/**
 * Quick helper to extract elements from POSCAR header (VASP 5 format line 6)
 */
function getPoscarElements(content: string): string[] {
    const lines = content.split(/\r?\n/);
    if (lines.length < 6) return [];

    // VASP 5: Line 6 contains element symbols
    // VASP 4: Line 6 contains counts (Implicit elements). We can't validate species names against VASP 4 POSCAR easily.

    const line6 = lines[5].trim();
    // Check if line 6 is numeric (VASP 4) or strings (VASP 5)
    // "1 14" vs "H H_pv"
    const tokens = line6.split(/\s+/);
    const isVasp5 = isNaN(parseFloat(tokens[0]));

    if (isVasp5) {
        return tokens;
    }
    return []; // Cannot validate against VASP 4 implicit POSCAR
}
