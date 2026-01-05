import { Diagnostic, DiagnosticSeverity, Range } from 'vscode-languageserver-types';
import { IncarDocument, IncarStatement } from './parsing';
import { VASP_TAGS, TagDefinition } from '../../data/vasp-tags';
import { isNumber, isInteger, createDiagnostic } from '../../utils/util';

/**
 * Validate a parsed INCAR document and return diagnostics.
 *
 * @param doc - Parsed INCAR document.
 * @returns Diagnostics describing invalid tags or values.
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
            diagnostics.push(
                createDiagnostic(
                    stmt.tag.range,
                    `Unknown tag '${tagName}'. Check spelling or version compatibility.`,
                    DiagnosticSeverity.Warning
                )
            );
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
                diagnostics.push(
                    createDiagnostic(
                        stmt.values[0].range,
                        'Warning: Time step is extremely fast (Presto). Please consider slowing down to Andante (0.5 - 2.0 fs) to ensure the ions can keep the beat.',
                        DiagnosticSeverity.Warning
                    )
                );
            }
        }
    }

    return diagnostics;
}

function validateValueType(stmt: IncarStatement, def: TagDefinition, diagnostics: Diagnostic[]) {
    // If the expected type is 'string', we treat all tokens as one string value.
    if (def.type === 'string') {
        // Combine tokens for "SYSTEM = Fe oxide"
        const fullValue = stmt.values.map((v) => v.text).join(' ');

        if (def.options) {
            // Options are usually single keywords.
            if (stmt.values.length > 1) {
                // Warn if multiple words provided for an Option type tag
                // e.g. ALGO = Fast Damped (Invalid)
                diagnostics.push(
                    createDiagnostic(
                        Range.create(stmt.values[1].range.start, stmt.values[stmt.values.length - 1].range.end),
                        `Tag '${stmt.tag.text}' expects a single keyword option, but multiple were found.`,
                        DiagnosticSeverity.Warning
                    )
                );
            } else {
                const found = def.options.find((opt) => opt.toLowerCase() === fullValue.toLowerCase());
                if (!found) {
                    diagnostics.push(
                        createDiagnostic(
                            stmt.values[0].range,
                            `Invalid option '${fullValue}'. Allowed: ${def.options.join(', ')}`,
                            DiagnosticSeverity.Error
                        )
                    );
                }
            }
        }
        return;
    }

    // For non-strings (int, float, bool), expect single value usually
    if (def.type !== 'array' && stmt.values.length > 1) {
        diagnostics.push(
            createDiagnostic(
                Range.create(stmt.values[1].range.start, stmt.values[stmt.values.length - 1].range.end),
                `Tag '${stmt.tag.text}' expects a single value, but multiple were found.`,
                DiagnosticSeverity.Warning
            )
        );
    }

    // Check the first value (primary value)
    const valToken = stmt.values[0];
    const valText = valToken.text;

    switch (def.type) {
        case 'int':
            if (!isInteger(valText)) {
                diagnostics.push(
                    createDiagnostic(
                        valToken.range,
                        `Expected integer for '${stmt.tag.text}', found '${valText}'.`,
                        DiagnosticSeverity.Error
                    )
                );
            }
            break;
        case 'float':
            if (!isNumber(valText)) {
                diagnostics.push(
                    createDiagnostic(
                        valToken.range,
                        `Expected number for '${stmt.tag.text}', found '${valText}'.`,
                        DiagnosticSeverity.Error
                    )
                );
            }
            break;
        case 'bool': {
            // VASP booleans: .TRUE., .FALSE., T, F, True, False
            const lower = valText.toLowerCase();
            const isBool =
                lower === '.true.' ||
                lower === '.false.' ||
                lower === 't' ||
                lower === 'f' ||
                lower === 'true' ||
                lower === 'false';

            if (!isBool) {
                diagnostics.push(
                    createDiagnostic(
                        valToken.range,
                        `Expected boolean (T/F/.TRUE./.FALSE.) for '${stmt.tag.text}', found '${valText}'.`,
                        DiagnosticSeverity.Error
                    )
                );
            }
            break;
        }
        case 'array':
            // Arrays: Check all items are valid numbers?
            // Usually arrays in VASP are numbers.
            for (const v of stmt.values) {
                if (!isNumber(v.text)) {
                    diagnostics.push(
                        createDiagnostic(
                            v.range,
                            `Expected number in array for '${stmt.tag.text}', found '${v.text}'.`,
                            DiagnosticSeverity.Error
                        )
                    );
                }
            }
            break;
    }
}
