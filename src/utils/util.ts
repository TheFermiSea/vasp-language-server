import { Diagnostic, DiagnosticSeverity, Range } from 'vscode-languageserver-types';

/**
 * Counts the elements in an array until a condition is met.
 * Useful for finding the 'valid' part of a line before invalid data starts.
 *
 * @param arr - The array to iterate over.
 * @param condition - The predicate. The count stops *before* the first element that satisfies this.
 * @returns The number of elements before the condition became true.
 */
export function countUntil<T>(arr: Array<T>, condition: (val: T) => boolean): number {
    let count = 0;
    while (count < arr.length && !condition(arr[count])) {
        count++;
    }
    return count;
}

/**
 * Build a Diagnostic object with a consistent default source/severity.
 *
 * @param range - Range where the diagnostic applies.
 * @param message - Human-readable message for the user.
 * @param severity - Severity level (default: Error).
 * @param code - Optional diagnostic code identifier.
 * @param source - Optional diagnostic source (default: vasp-lsp).
 * @returns A populated Diagnostic object.
 */
export function createDiagnostic(
    range: Range,
    message: string,
    severity: DiagnosticSeverity = DiagnosticSeverity.Error,
    code?: string,
    source: string = 'vasp-lsp'
): Diagnostic {
    const diagnostic: Diagnostic = {
        range,
        message,
        severity,
        source
    };

    if (code !== undefined) {
        diagnostic.code = code;
    }

    return diagnostic;
}

/**
 * Check whether a string matches a VASP-style numeric literal.
 *
 * Supports repetition prefixes (e.g. 5*1.2) and Fortran notation (e.g. 1.0D-3).
 *
 * @param str - Input string to validate.
 * @returns True when the string is a valid VASP numeric literal.
 */
export function isNumber(str: string): boolean {
    // Standard VASP number: optionally starts with "N*" repetition factor.
    // e.g., 5*1.2, 10*0, -1.2e-5, .5
    return /^(\d+\*)?[-+]?(\d+(\.\d*)?|\.\d+)([eEdD][-+]?\d+)?$/.test(str);
}

/**
 * Checks if a string represents a valid integer.
 *
 * @param str - The string to check.
 * @returns True if valid integer.
 */
export function isInteger(str: string): boolean {
    return /^-?\d+$/.test(str);
}

/**
 * Checks if a string consists only of alphabetic characters.
 * Used for validating species names (e.g., 'Fe', 'Li').
 *
 * @param str - The string to check.
 * @returns True if only letters.
 */
export function isLetters(str: string): boolean {
    return /^[a-zA-Z]+$/.test(str);
}

/**
 * Computes the Levenshtein distance between two strings.
 *
 * @param a - First string.
 * @param b - Second string.
 * @returns The edit distance between the two strings.
 */
export function levenshteinDistance(a: string, b: string): number {
    const matrix: number[][] = [];

    // Increment along the first column of each row
    for (let i = 0; i <= b.length; i++) {
        matrix[i] = [i];
    }

    // Increment each column in the first row
    for (let j = 0; j <= a.length; j++) {
        matrix[0][j] = j;
    }

    // Fill in the rest of the matrix
    for (let i = 1; i <= b.length; i++) {
        for (let j = 1; j <= a.length; j++) {
            if (b.charAt(i - 1) === a.charAt(j - 1)) {
                matrix[i][j] = matrix[i - 1][j - 1];
            } else {
                matrix[i][j] = Math.min(
                    matrix[i - 1][j - 1] + 1, // substitution
                    Math.min(
                        matrix[i][j - 1] + 1, // insertion
                        matrix[i - 1][j] + 1 // deletion
                    )
                );
            }
        }
    }

    return matrix[b.length][a.length];
}
