/**
 * CRYSTAL23 input file validation
 */

import { Diagnostic, DiagnosticSeverity } from 'vscode-languageserver-types';
import { CrystalDocument, CrystalStatement } from './parsing';
import { CRYSTAL_TAGS, CrystalTagDefinition } from '../../data/crystal-tags';
import { levenshteinDistance } from '../../utils/util';

/**
 * Validate a parsed CRYSTAL23 document
 */
export function validateCrystal(document: CrystalDocument): Diagnostic[] {
    const diagnostics: Diagnostic[] = [];

    // Add parse errors as diagnostics
    for (const error of document.errors) {
        diagnostics.push({
            range: error.range,
            severity: error.severity === 'error' ? DiagnosticSeverity.Error : DiagnosticSeverity.Warning,
            message: error.message,
            source: 'crystal-lsp'
        });
    }

    // Validate geometry type
    if (document.geometry.type === 'UNKNOWN') {
        diagnostics.push({
            range: document.geometry.typeToken?.range || {
                start: { line: 1, character: 0 },
                end: { line: 1, character: 10 }
            },
            severity: DiagnosticSeverity.Error,
            message: 'Invalid geometry type. Expected: CRYSTAL, SLAB, POLYMER, MOLECULE, HELIX, or EXTERNAL',
            source: 'crystal-lsp'
        });
    }

    // Validate space group for CRYSTAL type
    if (document.geometry.type === 'CRYSTAL' && document.geometry.spaceGroup !== undefined) {
        if (document.geometry.spaceGroup < 1 || document.geometry.spaceGroup > 230) {
            diagnostics.push({
                range: document.geometry.spaceGroupToken?.range || {
                    start: { line: 2, character: 0 },
                    end: { line: 2, character: 5 }
                },
                severity: DiagnosticSeverity.Error,
                message: `Invalid space group: ${document.geometry.spaceGroup}. Must be between 1 and 230`,
                source: 'crystal-lsp'
            });
        }
    }

    // Validate layer group for SLAB type
    if (document.geometry.type === 'SLAB' && document.geometry.spaceGroup !== undefined) {
        if (document.geometry.spaceGroup < 1 || document.geometry.spaceGroup > 80) {
            diagnostics.push({
                range: document.geometry.spaceGroupToken?.range || {
                    start: { line: 2, character: 0 },
                    end: { line: 2, character: 5 }
                },
                severity: DiagnosticSeverity.Warning,
                message: `Layer group ${document.geometry.spaceGroup} is outside typical range (1-80)`,
                source: 'crystal-lsp'
            });
        }
    }

    // Validate all statements
    for (const statement of document.allStatements) {
        const keyword = statement.keyword.text.toUpperCase();
        const tagDef = CRYSTAL_TAGS[keyword];

        if (!tagDef) {
            // Unknown keyword - suggest similar ones
            const suggestion = findSimilarKeyword(keyword);
            diagnostics.push({
                range: statement.keyword.range,
                severity: DiagnosticSeverity.Warning,
                message: suggestion
                    ? `Unknown keyword: ${keyword}. Did you mean '${suggestion}'?`
                    : `Unknown keyword: ${keyword}`,
                source: 'crystal-lsp',
                data: suggestion ? { suggestion } : undefined
            });
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
                diagnostics.push({
                    range: statement.range,
                    severity: DiagnosticSeverity.Error,
                    message: `${statement.keyword.text} expects ${min}-${max} arguments, got ${statement.values.length}`,
                    source: 'crystal-lsp'
                });
            }
        } else if (statement.values.length !== tagDef.argCount) {
            diagnostics.push({
                range: statement.range,
                severity: DiagnosticSeverity.Error,
                message: `${statement.keyword.text} expects ${tagDef.argCount} argument(s), got ${statement.values.length}`,
                source: 'crystal-lsp'
            });
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
            diagnostics.push({
                range: statement.values[i].range,
                severity: DiagnosticSeverity.Error,
                message: `Expected integer for argument ${i + 1}, got '${value}'`,
                source: 'crystal-lsp'
            });
        } else if (expectedType === 'float' && !/^-?\d*\.?\d+([eEdD][+-]?\d+)?$/.test(value)) {
            diagnostics.push({
                range: statement.values[i].range,
                severity: DiagnosticSeverity.Error,
                message: `Expected number for argument ${i + 1}, got '${value}'`,
                source: 'crystal-lsp'
            });
        }
    }
}

function validateKeywordSpecific(
    statement: CrystalStatement,
    tagDef: CrystalTagDefinition,
    diagnostics: Diagnostic[]
): void {
    const keyword = statement.keyword.text.toUpperCase();

    // SHRINK validation - ensure reasonable grid sizes
    if (keyword === 'SHRINK' && statement.values.length >= 2) {
        const is = parseInt(statement.values[0].text, 10);
        const ip = parseInt(statement.values[1].text, 10);

        if (!isNaN(is) && is < 1) {
            diagnostics.push({
                range: statement.values[0].range,
                severity: DiagnosticSeverity.Error,
                message: `SHRINK IS must be >= 1, got ${is}`,
                source: 'crystal-lsp'
            });
        }
        if (!isNaN(ip) && ip < 1) {
            diagnostics.push({
                range: statement.values[1].range,
                severity: DiagnosticSeverity.Error,
                message: `SHRINK IP (Gilat net) must be >= 1, got ${ip}`,
                source: 'crystal-lsp'
            });
        }
        if (!isNaN(is) && !isNaN(ip) && ip > is) {
            diagnostics.push({
                range: statement.range,
                severity: DiagnosticSeverity.Warning,
                message: `SHRINK IP (${ip}) is typically <= IS (${is})`,
                source: 'crystal-lsp'
            });
        }
    }

    // TOLINTEG validation - 5 integers required
    if (keyword === 'TOLINTEG' && statement.values.length > 0) {
        if (statement.values.length !== 5) {
            diagnostics.push({
                range: statement.range,
                severity: DiagnosticSeverity.Error,
                message: `TOLINTEG requires exactly 5 integers (overlap/penetration thresholds)`,
                source: 'crystal-lsp'
            });
        }
    }

    // FMIXING validation - percentage 0-100
    if (keyword === 'FMIXING' && statement.values.length >= 1) {
        const pct = parseInt(statement.values[0].text, 10);
        if (!isNaN(pct) && (pct < 0 || pct > 100)) {
            diagnostics.push({
                range: statement.values[0].range,
                severity: DiagnosticSeverity.Warning,
                message: `FMIXING percentage should be 0-100, got ${pct}`,
                source: 'crystal-lsp'
            });
        }
    }

    // MAXCYCLE validation
    if (keyword === 'MAXCYCLE' && statement.values.length >= 1) {
        const cycles = parseInt(statement.values[0].text, 10);
        if (!isNaN(cycles) && cycles < 1) {
            diagnostics.push({
                range: statement.values[0].range,
                severity: DiagnosticSeverity.Error,
                message: `MAXCYCLE must be >= 1`,
                source: 'crystal-lsp'
            });
        }
        if (!isNaN(cycles) && cycles > 500) {
            diagnostics.push({
                range: statement.values[0].range,
                severity: DiagnosticSeverity.Warning,
                message: `MAXCYCLE=${cycles} is unusually high. SCF issues may need different approach`,
                source: 'crystal-lsp'
            });
        }
    }
}

function validateRequiredKeywords(document: CrystalDocument, diagnostics: Diagnostic[]): void {
    const keywords = new Set(document.allStatements.map((s) => s.keyword.text.toUpperCase()));

    // SHRINK is almost always required
    if (!keywords.has('SHRINK') && document.geometry.type !== 'MOLECULE') {
        diagnostics.push({
            range: { start: { line: 0, character: 0 }, end: { line: 0, character: 1 } },
            severity: DiagnosticSeverity.Warning,
            message: 'SHRINK keyword not found. K-point sampling is typically required for periodic systems',
            source: 'crystal-lsp'
        });
    }

    // Check END keywords balance
    const endCount = document.allTokens.filter((t) => t.type === 'end').length;
    if (endCount < 1) {
        diagnostics.push({
            range: {
                start: { line: document.lines.length - 1, character: 0 },
                end: { line: document.lines.length - 1, character: 1 }
            },
            severity: DiagnosticSeverity.Warning,
            message: 'Missing END keyword. CRYSTAL23 input typically ends with END',
            source: 'crystal-lsp'
        });
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
