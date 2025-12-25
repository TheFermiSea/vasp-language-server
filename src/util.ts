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
 * Checks if a string represents a valid real number (float or integer).
 * Matches scientific notation (e.g. 1.0E-5).
 * 
 * @param str - The string to check.
 * @returns True if valid number.
 */
export function isNumber(str: string): boolean {
    return /^-?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?$/.test(str);
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
