import { Diagnostic, DiagnosticSeverity, Range } from "vscode-languageserver-types";
import { IncarDocument, IncarStatement } from "./incar-parsing";
import { VASP_TAGS, TagDefinition } from "./data/vasp-tags";
import { isNumber, isInteger } from "./util";

/**
 * Validates a parsed INCAR document.
 */
export function validateIncar(doc: IncarDocument): Diagnostic[] {
    const diagnostics: Diagnostic[] = [];

    for (const stmt of doc.statements) {
        // 1. Add any parsing errors (e.g. missing equals)
        diagnostics.push(...stmt.parsingErrors);

        const tagName = stmt.tag.text.toUpperCase();
        const definition = VASP_TAGS[tagName];

        // 2. Unknown Tag Check
        if (!definition) {
            diagnostics.push({
                message: `Unknown tag '${tagName}'. Check spelling or version compatibility.`,
                range: stmt.tag.range,
                severity: DiagnosticSeverity.Warning,
                source: "VASP"
            });
            continue;
        }

        // 3. Type Validation
        if (stmt.values.length > 0) {
            validateValueType(stmt, definition, diagnostics);
        }

        // --- Easter Egg: POTIM Tempo Check ---
        if (tagName === 'POTIM') {
            const val = parseFloat(stmt.values[0].text);
            if (!isNaN(val) && val > 3.0) {
                diagnostics.push({
                    message: `Warning: Time step is extremely fast (Presto). Please consider slowing down to Andante (0.5 - 2.0 fs) to ensure the ions can keep the beat.`,
                    range: stmt.values[0].range,
                    severity: DiagnosticSeverity.Warning,
                    source: "Maestro"
                });
            }
        }
    }

    return diagnostics;
}

function validateValueType(stmt: IncarStatement, def: TagDefinition, diagnostics: Diagnostic[]) {
    // If the expected type is 'string', we treat all tokens as one string value.
    if (def.type === 'string') {
        // Combine tokens for "SYSTEM = Fe oxide"
        const fullValue = stmt.values.map(v => v.text).join(" ");

        if (def.options) {
            // Options are usually single keywords.
            if (stmt.values.length > 1) {
                // Warn if multiple words provided for an Option type tag
                // e.g. ALGO = Fast Damped (Invalid)
                diagnostics.push({
                    message: `Tag '${stmt.tag.text}' expects a single keyword option, but multiple were found.`,
                    range: Range.create(stmt.values[1].range.start, stmt.values[stmt.values.length - 1].range.end),
                    severity: DiagnosticSeverity.Warning,
                    source: "VASP"
                });
            } else {
                const found = def.options.find(opt => opt.toLowerCase() === fullValue.toLowerCase());
                if (!found) {
                    diagnostics.push({
                        message: `Invalid option '${fullValue}'. Allowed: ${def.options.join(', ')}`,
                        range: stmt.values[0].range,
                        severity: DiagnosticSeverity.Error,
                        source: "VASP"
                    });
                }
            }
        }
        return;
    }

    // For non-strings (int, float, bool), expect single value usually
    if (def.type !== 'array' && stmt.values.length > 1) {
        diagnostics.push({
            message: `Tag '${stmt.tag.text}' expects a single value, but multiple were found.`,
            range: Range.create(stmt.values[1].range.start, stmt.values[stmt.values.length - 1].range.end),
            severity: DiagnosticSeverity.Warning,
            source: "VASP"
        });
    }

    // Check the first value (primary value)
    const valToken = stmt.values[0];
    const valText = valToken.text;

    switch (def.type) {
        case 'int':
            if (!isInteger(valText)) {
                diagnostics.push({
                    message: `Expected integer for '${stmt.tag.text}', found '${valText}'.`,
                    range: valToken.range,
                    severity: DiagnosticSeverity.Error,
                    source: "VASP"
                });
            }
            break;
        case 'float':
            if (!isNumber(valText)) {
                diagnostics.push({
                    message: `Expected number for '${stmt.tag.text}', found '${valText}'.`,
                    range: valToken.range,
                    severity: DiagnosticSeverity.Error,
                    source: "VASP"
                });
            }
            break;
        case 'bool':
            // VASP booleans: .TRUE., .FALSE., T, F, True, False
            const lower = valText.toLowerCase();
            const isBool =
                lower === '.true.' || lower === '.false.' ||
                lower === 't' || lower === 'f' ||
                lower === 'true' || lower === 'false';

            if (!isBool) {
                diagnostics.push({
                    message: `Expected boolean (T/F/.TRUE./.FALSE.) for '${stmt.tag.text}', found '${valText}'.`,
                    range: valToken.range,
                    severity: DiagnosticSeverity.Error,
                    source: "VASP"
                });
            }
            break;
        case 'array':
            // Arrays: Check all items are valid numbers?
            // Usually arrays in VASP are numbers.
            for (const v of stmt.values) {
                if (!isNumber(v.text)) {
                    diagnostics.push({
                        message: `Expected number in array for '${stmt.tag.text}', found '${v.text}'.`,
                        range: v.range,
                        severity: DiagnosticSeverity.Error,
                        source: "VASP"
                    });
                }
            }
            break;
    }
}
