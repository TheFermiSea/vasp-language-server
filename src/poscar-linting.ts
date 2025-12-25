import { Diagnostic, DiagnosticSeverity, Range } from 'vscode-languageserver-types';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { isNumber, isInteger, isLetters } from './util';
import { parsePoscar, PoscarBlockType, PoscarLine, PoscarDocument } from './poscar-parsing';
import { countUntil } from './util';

/**
 * Main validation function for POSCAR files.
 * Orchestrates the parsing and applies specific linting rules to each block.
 *
 * @param document - The text document to validate.
 * @param parsed - Optional pre-parsed document.
 * @returns Array of Diagnostics (errors/warnings).
 */
export function validatePoscar(document: TextDocument, parsed?: PoscarDocument): Diagnostic[] {
    const poscarDoc = parsed || parsePoscar(document);
    const poscarLines = poscarDoc.lines;

    // 1. Block-level validation
    // Apply specific linter rules to each line based on its parsed type.
    const diagnostics = poscarLines.flatMap((l) => {
        const linter = (poscarBlockLinters as any)[l.type];
        return linter ? linter(l) : [];
    });

    // 2. Cross-block validation
    // Example: Check if atomic species names match the number of atom counts.
    // (VASP requires one count per species)
    const speciesNamesLine = poscarLines.find((l) => l.type === 'speciesNames');
    const numAtomsLine = poscarLines.find((l) => l.type === 'numAtoms');

    if (speciesNamesLine && numAtomsLine) {
        // countUntilComment returns the number of valid tokens before a comment start
        if (countUntilComment(speciesNamesLine) !== countUntilComment(numAtomsLine)) {
            diagnostics.push(
                createDiagnostic(
                    'Number of atoms must be specified for each atomic species.',
                    numAtomsLine.line.range,
                    DiagnosticSeverity.Error
                )
            );
        }
    }

    return diagnostics;
}

/**
 * Signature for a linter function that checks a specific line type.
 */
type Linter = (poscarLine: PoscarLine) => Diagnostic[];

/**
 * Map of BlockTypes to their specific Linter functions.
 */
const poscarBlockLinters: Readonly<Record<PoscarBlockType, Linter>> = {
    // Comments are free-form, no linting.
    comment: () => [],

    scaling: (poscarLine) => {
        const diagnostics: Diagnostic[] = [];
        if (isEmptyLine(poscarLine, diagnostics)) {
            return diagnostics;
        }

        const tokens = poscarLine.tokens;
        const numVals = countUntilComment(poscarLine, diagnostics);

        // Scaling factor must be 1 (universal) or 3 (per-axis).
        if (numVals === 2 || numVals > 3) {
            const startToken = tokens[0];
            const endToken = tokens[numVals - 1];
            // Highlight the whole invalid range
            const unionRange = Range.create(startToken.range.start, endToken.range.end);

            diagnostics.push(
                createDiagnostic(
                    'The number of scaling factors must be either 1 or 3.',
                    unionRange,
                    DiagnosticSeverity.Error
                )
            );
        } else if (numVals === 3) {
            // Per-axis scaling factors must be positive.
            tokens.forEach((t) => {
                if (+t.text < 0) {
                    diagnostics.push(
                        createDiagnostic(
                            'Individual scaling factors must be positive.',
                            t.range,
                            DiagnosticSeverity.Error
                        )
                    );
                }
            });
        }
        return diagnostics;
    },

    lattice: lintVector,

    speciesNames: (poscarLine) => {
        // Validate that all species names are just letters
        const diagnostics = poscarLine.tokens
            .filter((t: any) => !isLetters(t.text))
            .map((t: any) => createDiagnostic(`Species name '${t.text}' is invalid.`, t.range, DiagnosticSeverity.Error));
        isEmptyLine(poscarLine, diagnostics);
        return diagnostics;
    },

    numAtoms: (poscarLine) => {
        const diagnostics: Diagnostic[] = [];
        if (isEmptyLine(poscarLine, diagnostics)) {
            return diagnostics;
        }
        const numVals = countUntilComment(poscarLine, diagnostics);

        // Atom counts must be positive integers
        poscarLine.tokens
            .slice(0, numVals)
            .filter((t) => !isInteger(t.text) || +t.text <= 0)
            .forEach((t) => {
                diagnostics.push(
                    createDiagnostic(
                        'Number of atoms needs to be a positive integer.',
                        t.range,
                        DiagnosticSeverity.Error
                    )
                );
            });
        return diagnostics;
    },

    selDynamics: (poscarLine) => lintConstLine(poscarLine, 'selective dynamics'),

    positionMode: (poscarLine) => lintMode(poscarLine, 'direct'),

    positions: lintVector,

    positionsSelDyn: (poscarLine) => {
        // First lint the coordinates (first 3 numbers)
        const diagnostics = lintVector(poscarLine);
        const tokens = poscarLine.tokens;

        // Then lint the flags (indices 3, 4, 5)
        tokens
            .slice(3, 6)
            .filter((t) => t.text !== 'T' && t.text !== 'F')
            .forEach((t) => {
                diagnostics.push(
                    createDiagnostic(
                        "Selective dynamics flag must be either 'T' or 'F'.",
                        t.range,
                        DiagnosticSeverity.Error
                    )
                );
            });

        // Ensure we have exactly 3 flags
        if (tokens.length <= 3) {
            diagnostics.push(
                createDiagnostic(
                    'There must be 3 selective-dynamics flags. Too few given.',
                    Range.create(poscarLine.line.range.end, poscarLine.line.range.end),
                    DiagnosticSeverity.Error
                )
            );
        } else if (tokens.length < 6) {
            const startToken = tokens[3];
            const endToken = tokens[tokens.length - 1];
            const unionRange = Range.create(startToken.range.start, endToken.range.end);
            diagnostics.push(
                createDiagnostic(
                    'There must be 3 selective-dynamics flags. Too few given.',
                    unionRange,
                    DiagnosticSeverity.Error
                )
            );
        }
        return diagnostics;
    },

    lattVelocitiesStart: (poscarLine) => lintConstLine(poscarLine, 'Lattice velocities and vectors'),

    lattVelocitiesState: (poscarLine) => {
        const tokens = poscarLine.tokens;
        const diagnostics: Diagnostic[] = [];
        if (!isEmptyLine(poscarLine, diagnostics) && !isInteger(tokens[0].text)) {
            diagnostics.push(
                createDiagnostic(
                    'Initialization state needs to be an integer',
                    tokens[0].range,
                    DiagnosticSeverity.Error
                )
            );
        }
        return diagnostics;
    },

    lattVelocitiesVels: lintVector,
    lattVelocitiesLatt: lintVector,
    velocityMode: (poscarLine) => lintMode(poscarLine),
    velocities: lintVector
};

/**
 * Helper to create a standardized Diagnostic object.
 */
function createDiagnostic(message: string, range: Range, severity: DiagnosticSeverity): Diagnostic {
    return {
        message: message,
        range: range,
        severity: severity,
        source: 'VASP support'
    };
}

/**
 * Checks if a line is completely empty (no tokens).
 * If so, optionally adds an Error diagnostic.
 */
function isEmptyLine(poscarLine: PoscarLine, diagnostics?: Diagnostic[]): boolean {
    if (poscarLine.tokens.length === 0) {
        diagnostics?.push(
            createDiagnostic(
                'Line must not be empty.',
                poscarLine.line.rangeIncludingLineBreak,
                DiagnosticSeverity.Error
            )
        );
        return true;
    }
    return false;
}

/**
 * Counts valid tokens until a comment starts.
 * Also lints for implicit/ambiguous comments (garbage without #).
 */
function countUntilComment(poscarLine: PoscarLine, diagnostics?: Diagnostic[]): number {
    const tokens = poscarLine.tokens;
    const numVals = countUntil(tokens, (t: any) => t.type === 'comment');

    // If there is "trailing garbage" that isn't explicitly marked as a comment (no # or !), warn the user.
    if (diagnostics && numVals < tokens.length && !/^[#!]/.test(tokens[numVals].text)) {
        const tokenRange = tokens[numVals].range;
        const endOfLine = poscarLine.line.range.end;
        const range = Range.create(tokenRange.start, endOfLine);

        diagnostics.push(
            createDiagnostic(
                'The remainder of this line is ignored by VASP. ' +
                "Consider placing a '#' or '!' in front to make the intention clearer.",
                range,
                DiagnosticSeverity.Warning
            )
        );
    }
    return numVals;
}

/**
 * Standard linter for vector lines (e.g., lattice vectors, positions).
 * Expects exactly 3 numbers.
 */
function lintVector(poscarLine: PoscarLine): Diagnostic[] {
    const diagnostics: Diagnostic[] = [];
    if (isEmptyLine(poscarLine, diagnostics)) {
        return diagnostics;
    }

    const tokens = poscarLine.tokens;
    countUntilComment(poscarLine, diagnostics);

    // Check first 3 tokens are numbers
    tokens.slice(0, 3).forEach((t) => {
        if (!isNumber(t.text)) {
            diagnostics.push(createDiagnostic('Vector component must be a number.', t.range, DiagnosticSeverity.Error));
        }
    });

    // Ensure we have at least 3
    if (tokens.length < 3) {
        const startToken = tokens[0];
        const endToken = tokens[tokens.length - 1];
        const unionRange = Range.create(startToken.range.start, endToken.range.end);

        diagnostics.push(
            createDiagnostic('Vector must consist of 3 numbers. Too few given.', unionRange, DiagnosticSeverity.Error)
        );
    }
    return diagnostics;
}

/**
 * Linter for constant/keyword lines.
 * Checks if the keyword starts with the expected character.
 */
function lintConstLine(poscarLine: PoscarLine, content: string): Diagnostic[] {
    const diagnostics: Diagnostic[] = [];
    if (isEmptyLine(poscarLine, diagnostics)) {
        return diagnostics;
    }

    const token = poscarLine.tokens[0];
    const startLower = content[0].toLowerCase();
    const startUpper = startLower.toUpperCase();
    const regex = new RegExp(`^[${startLower}${startUpper}]`);

    if (!regex.test(token.text)) {
        diagnostics.push(
            createDiagnostic(
                `First non-space character on line must be '${startLower}' or '${startUpper}'.`,
                token.range,
                DiagnosticSeverity.Error
            )
        );
    } else if (!content.startsWith(token.text.toLowerCase())) {
        // Warn if it matches the letter but isn't the full expected word
        diagnostics.push(
            createDiagnostic(
                `Consider specifying '${content}' to avoid potential mistakes.`,
                token.range,
                DiagnosticSeverity.Warning
            )
        );
    }
    return diagnostics;
}

/**
 * Special linter for Mode lines (Direct vs Cartesian).
 */
function lintMode(poscarLine: PoscarLine, emptyMode?: string) {
    const diagnostics: Diagnostic[] = [];
    if (poscarLine.tokens.length === 0) {
        if (emptyMode) {
            diagnostics.push(
                createDiagnostic(
                    `Consider specifying '${emptyMode}' instead of an empty line to avoid potential mistakes.`,
                    poscarLine.line.rangeIncludingLineBreak,
                    DiagnosticSeverity.Warning
                )
            );
        }
    } else {
        const token = poscarLine.tokens[0];
        const firstLetter = token.text ? token.text[0].toLowerCase() : '';

        // VASP allows 'C', 'K', or 'D' to specify modes.
        if (firstLetter === 'c' || firstLetter === 'k') {
            if (!'artesian'.startsWith(token.text.slice(1).toLowerCase())) {
                diagnostics.push(
                    createDiagnostic(
                        "Consider specifying 'cartesian' to avoid potential mistakes.",
                        token.range,
                        DiagnosticSeverity.Warning
                    )
                );
            }
        } else {
            if (!'direct'.startsWith(token.text.toLowerCase())) {
                diagnostics.push(
                    createDiagnostic(
                        "Consider specifying 'direct' to avoid potential mistakes.",
                        token.range,
                        DiagnosticSeverity.Warning
                    )
                );
            }
        }
    }
    return diagnostics;
}
