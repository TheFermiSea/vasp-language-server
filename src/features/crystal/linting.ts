/**
 * CRYSTAL23 input file validation
 */

import { Diagnostic, DiagnosticSeverity, Range } from 'vscode-languageserver-types';
import { CrystalDocument, CrystalStatement } from './parsing';
import { CRYSTAL_TAGS, CrystalTagDefinition } from '../../data/crystal-tags';
import { createDiagnostic, levenshteinDistance } from '../../utils/util';

/**
 * Validate a parsed CRYSTAL23 document.
 *
 * @param document - Parsed CRYSTAL document.
 * @returns Diagnostics for invalid geometry, keywords, or arguments.
 */
export function validateCrystal(document: CrystalDocument): Diagnostic[] {
    const diagnostics: Diagnostic[] = [];

    // Add parse errors as diagnostics
    for (const error of document.errors) {
        diagnostics.push(
            createCrystalDiagnostic(
                error.range,
                error.message,
                error.severity === 'error' ? DiagnosticSeverity.Error : DiagnosticSeverity.Warning
            )
        );
    }

    // Validate geometry type
    if (document.geometry.type === 'UNKNOWN') {
        diagnostics.push(
            createCrystalDiagnostic(
                document.geometry.typeToken?.range || {
                    start: { line: 1, character: 0 },
                    end: { line: 1, character: 10 }
                },
                'Invalid geometry type. Expected: CRYSTAL, SLAB, POLYMER, MOLECULE, HELIX, or EXTERNAL',
                DiagnosticSeverity.Error
            )
        );
    }

    // Validate space group for CRYSTAL type
    if (document.geometry.type === 'CRYSTAL' && document.geometry.spaceGroup !== undefined) {
        if (document.geometry.spaceGroup < 1 || document.geometry.spaceGroup > 230) {
            diagnostics.push(
                createCrystalDiagnostic(
                    document.geometry.spaceGroupToken?.range || {
                        start: { line: 2, character: 0 },
                        end: { line: 2, character: 5 }
                    },
                    `Invalid space group: ${document.geometry.spaceGroup}. Must be between 1 and 230`,
                    DiagnosticSeverity.Error
                )
            );
        }
    }

    // Validate layer group for SLAB type
    if (document.geometry.type === 'SLAB' && document.geometry.spaceGroup !== undefined) {
        if (document.geometry.spaceGroup < 1 || document.geometry.spaceGroup > 80) {
            diagnostics.push(
                createCrystalDiagnostic(
                    document.geometry.spaceGroupToken?.range || {
                        start: { line: 2, character: 0 },
                        end: { line: 2, character: 5 }
                    },
                    `Layer group ${document.geometry.spaceGroup} is outside typical range (1-80)`,
                    DiagnosticSeverity.Warning
                )
            );
        }
    }

    // Validate all statements
    for (const statement of document.allStatements) {
        const keyword = statement.keyword.text.toUpperCase();
        const tagDef = CRYSTAL_TAGS[keyword];

        if (!tagDef) {
            // Unknown keyword - suggest similar ones
            const suggestion = findSimilarKeyword(keyword);
            const diagnostic = createCrystalDiagnostic(
                statement.keyword.range,
                suggestion
                    ? `Unknown keyword: ${keyword}. Did you mean '${suggestion}'?`
                    : `Unknown keyword: ${keyword}`,
                DiagnosticSeverity.Warning
            );
            if (suggestion) {
                diagnostic.data = { suggestion };
            }
            diagnostics.push(diagnostic);
            continue;
        }

        // Validate argument count
        validateArgumentCount(statement, tagDef, diagnostics);

        // Validate argument types
        validateArgumentTypes(statement, tagDef, diagnostics);

        // Keyword-specific validation
        validateKeywordSpecific(statement, tagDef, diagnostics);
    }

    // Check for required keywords
    validateRequiredKeywords(document, diagnostics);

    return diagnostics;
}

function createCrystalDiagnostic(
    range: Range,
    message: string,
    severity: DiagnosticSeverity,
    code?: string
): Diagnostic {
    return createDiagnostic(range, message, severity, code, 'crystal-lsp');
}

function validateArgumentCount(
    statement: CrystalStatement,
    tagDef: CrystalTagDefinition,
    diagnostics: Diagnostic[]
): void {
    if (tagDef.argCount !== undefined) {
        if (Array.isArray(tagDef.argCount)) {
            // Variable argument count
            const [min, max] = tagDef.argCount;
            if (statement.values.length < min || statement.values.length > max) {
                diagnostics.push(
                    createCrystalDiagnostic(
                        statement.range,
                        `${statement.keyword.text} expects ${min}-${max} arguments, got ${statement.values.length}`,
                        DiagnosticSeverity.Error
                    )
                );
            }
        } else if (statement.values.length !== tagDef.argCount) {
            diagnostics.push(
                createCrystalDiagnostic(
                    statement.range,
                    `${statement.keyword.text} expects ${tagDef.argCount} argument(s), got ${statement.values.length}`,
                    DiagnosticSeverity.Error
                )
            );
        }
    }
}

function validateArgumentTypes(
    statement: CrystalStatement,
    tagDef: CrystalTagDefinition,
    diagnostics: Diagnostic[]
): void {
    if (!tagDef.argTypes) return;

    for (let i = 0; i < statement.values.length && i < tagDef.argTypes.length; i++) {
        const expectedType = tagDef.argTypes[i];
        const value = statement.values[i].text;

        if (expectedType === 'int' && !/^-?\d+$/.test(value)) {
            diagnostics.push(
                createCrystalDiagnostic(
                    statement.values[i].range,
                    `Expected integer for argument ${i + 1}, got '${value}'`,
                    DiagnosticSeverity.Error
                )
            );
        } else if (expectedType === 'float' && !/^-?\d*\.?\d+([eEdD][+-]?\d+)?$/.test(value)) {
            diagnostics.push(
                createCrystalDiagnostic(
                    statement.values[i].range,
                    `Expected number for argument ${i + 1}, got '${value}'`,
                    DiagnosticSeverity.Error
                )
            );
        }
    }
}

function validateKeywordSpecific(
    statement: CrystalStatement,
    tagDef: CrystalTagDefinition,
    diagnostics: Diagnostic[]
): void {
    void tagDef;

    const keyword = statement.keyword.text.toUpperCase();

    // SHRINK validation - ensure reasonable grid sizes
    if (keyword === 'SHRINK' && statement.values.length >= 2) {
        const is = parseInt(statement.values[0].text, 10);
        const ip = parseInt(statement.values[1].text, 10);

        if (!isNaN(is) && is < 1) {
            diagnostics.push(
                createCrystalDiagnostic(
                    statement.values[0].range,
                    `SHRINK IS must be >= 1, got ${is}`,
                    DiagnosticSeverity.Error
                )
            );
        }
        if (!isNaN(ip) && ip < 1) {
            diagnostics.push(
                createCrystalDiagnostic(
                    statement.values[1].range,
                    `SHRINK IP (Gilat net) must be >= 1, got ${ip}`,
                    DiagnosticSeverity.Error
                )
            );
        }
        if (!isNaN(is) && !isNaN(ip) && ip > is) {
            diagnostics.push(
                createCrystalDiagnostic(
                    statement.range,
                    `SHRINK IP (${ip}) is typically <= IS (${is})`,
                    DiagnosticSeverity.Warning
                )
            );
        }
    }

    // TOLINTEG validation - 5 integers required
    if (keyword === 'TOLINTEG' && statement.values.length > 0) {
        if (statement.values.length !== 5) {
            diagnostics.push(
                createCrystalDiagnostic(
                    statement.range,
                    'TOLINTEG requires exactly 5 integers (overlap/penetration thresholds)',
                    DiagnosticSeverity.Error
                )
            );
        }
    }

    // FMIXING validation - percentage 0-100
    if (keyword === 'FMIXING' && statement.values.length >= 1) {
        const pct = parseInt(statement.values[0].text, 10);
        if (!isNaN(pct) && (pct < 0 || pct > 100)) {
            diagnostics.push(
                createCrystalDiagnostic(
                    statement.values[0].range,
                    `FMIXING percentage should be 0-100, got ${pct}`,
                    DiagnosticSeverity.Warning
                )
            );
        }
    }

    // MAXCYCLE validation
    if (keyword === 'MAXCYCLE' && statement.values.length >= 1) {
        const cycles = parseInt(statement.values[0].text, 10);
        if (!isNaN(cycles) && cycles < 1) {
            diagnostics.push(
                createCrystalDiagnostic(statement.values[0].range, 'MAXCYCLE must be >= 1', DiagnosticSeverity.Error)
            );
        }
        if (!isNaN(cycles) && cycles > 500) {
            diagnostics.push(
                createCrystalDiagnostic(
                    statement.values[0].range,
                    `MAXCYCLE=${cycles} is unusually high. SCF issues may need different approach`,
                    DiagnosticSeverity.Warning
                )
            );
        }
    }
}

function validateRequiredKeywords(document: CrystalDocument, diagnostics: Diagnostic[]): void {
    const keywords = new Set(document.allStatements.map((s) => s.keyword.text.toUpperCase()));

    // SHRINK is almost always required
    if (!keywords.has('SHRINK') && document.geometry.type !== 'MOLECULE') {
        diagnostics.push(
            createCrystalDiagnostic(
                { start: { line: 0, character: 0 }, end: { line: 0, character: 1 } },
                'SHRINK keyword not found. K-point sampling is typically required for periodic systems',
                DiagnosticSeverity.Warning
            )
        );
    }

    // Check END keywords balance
    const endCount = document.allTokens.filter((t) => t.type === 'end').length;
    if (endCount < 1) {
        diagnostics.push(
            createCrystalDiagnostic(
                {
                    start: { line: document.lines.length - 1, character: 0 },
                    end: { line: document.lines.length - 1, character: 1 }
                },
                'Missing END keyword. CRYSTAL23 input typically ends with END',
                DiagnosticSeverity.Warning
            )
        );
    }
}

function findSimilarKeyword(unknown: string): string | undefined {
    const allKeywords = Object.keys(CRYSTAL_TAGS);
    let bestMatch: string | undefined;
    let bestDistance = 4; // Threshold

    for (const keyword of allKeywords) {
        const dist = levenshteinDistance(unknown.toLowerCase(), keyword.toLowerCase());
        if (dist < bestDistance) {
            bestDistance = dist;
            bestMatch = keyword;
        }
    }

    return bestMatch;
}
