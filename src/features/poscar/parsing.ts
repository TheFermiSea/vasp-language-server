import { Range } from 'vscode-languageserver-types';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { countUntil, isNumber, isInteger, isLetters } from '../../utils/util';

/**
 * Enumeration of all possible token types in a POSCAR file.
 * Used for syntax highlighting (in future) and linting.
 */
export const tokenTypes = ['comment', 'string', 'number', 'constant', 'invalid'] as const;

type TokenType = (typeof tokenTypes)[number];

/**
 * Enumeration of the different logical blocks/lines found in a POSCAR file.
 * The order roughly corresponds to the VASP 5.x format specification.
 */
export type PoscarBlockType =
    | 'comment' // Comment line (first line)
    | 'scaling' // Scaling factor (Universal or per-axis)
    | 'lattice' // Lattice vectors (3 lines)
    | 'speciesNames' // Atomic species names (e.g., Fe O)
    | 'numAtoms' // Number of atoms per species
    | 'selDynamics' // 'Selective dynamics' switch (Optional)
    | 'positionMode' // 'Direct' or 'Cartesian'
    | 'positions' // Atomic positions (without flags)
    | 'positionsSelDyn' // Atomic positions (with T/F flags)
    | 'lattVelocitiesStart' // Start of Lattice Velocities (Optional)
    | 'lattVelocitiesState' // State integer for lattice velocities
    | 'lattVelocitiesVels' // Lattice velocities vectors
    | 'lattVelocitiesLatt' // Lattice vectors (again? context dependent)
    | 'velocityMode' // 'Direct' or 'Cartesian' for velocities
    | 'velocities'; // Atom velocities

/**
 * A simplified interface representing a line of text in the LSP environment.
 * Mirrors `vscode.TextLine` but decoupled from the VS Code API.
 */
export interface LSPTextLine {
    text: string;
    lineNumber: number;
    range: Range;
    rangeIncludingLineBreak: Range;
    firstNonWhitespaceCharacterIndex: number;
    isEmptyOrWhitespace: boolean;
}

/**
 * Represents a fully parsed line in a POSCAR file.
 * Contains the raw text, the structural type of the line, and the breakdown of tokens.
 */
export interface PoscarLine {
    type: PoscarBlockType;
    tokens: Token[];
    line: LSPTextLine;
}

/**
 * Represents the full parsed result of a POSCAR file.
 */
export interface PoscarDocument {
    lines: PoscarLine[];
}

/**
 * A single lexical token (e.g., a number, a keyword, a comment).
 */
interface Token {
    type?: TokenType;
    range: Range;
    text: string;
}

/**
 * Signature for a function that parses a raw text line into tokens.
 */
type Tokenizer = (line: LSPTextLine) => Token[];

/**
 * Map of BlockTypes to their specific Tokenizer functions.
 * Defines the grammar of each line type.
 */
const tokenizers: Readonly<Record<PoscarBlockType, Tokenizer>> = {
    // just tokenizes everything as a comment
    comment: (line) => {
        return tokenizeLine(line).map((t) => {
            t.type = 'comment';
            return t;
        });
    },
    // Scaling factor line. Can be 1 number (universal) or 3 numbers (axes).
    scaling: (line) => {
        const tokens = tokenizeLine(line);
        // Find first non-number to detect comments or invalid trailing garbage
        const numVals = countUntil(tokens, (t) => !isNumber(t.text));

        // Note: 3 separate scaling factors is a rare/legacy VASP feature
        if (numVals === 3) {
            tokens.slice(0, 3).forEach((t) => (t.type = +t.text > 0 ? 'number' : 'invalid'));
        } else {
            // Otherwise, expect 1 number. Anything else before comment is invalid.
            tokens.slice(0, numVals).forEach((t, tIdx) => {
                t.type = tIdx < 3 ? 'number' : 'invalid';
            });
        }
        // Mark remainder as comment
        tokens.slice(numVals).forEach((t) => (t.type = 'comment'));

        return tokens;
    },
    lattice: tokenizeVector,
    speciesNames: (line) => {
        return tokenizeLine(line).map((t) => {
            t.type = isLetters(t.text) ? 'string' : 'invalid';
            return t;
        });
    },
    numAtoms: (line) => {
        const tokens = tokenizeLine(line);
        const numVals = countUntil(tokens, (t) => !isInteger(t.text));

        tokens.slice(0, numVals).forEach((t) => (t.type = 'number'));
        tokens.slice(numVals).forEach((t) => (t.type = 'comment'));
        return tokens;
    },
    // Checks for 'Selective dynamics' (starts with 's' or 'S')
    selDynamics: (line) => tokenizeConstLine(line, /^[sS]/),
    positionMode: tokenizeConstLine,
    positions: tokenizeVector,
    // Positions with T/F flags for selective dynamics
    positionsSelDyn: (line) => {
        return tokenizeLine(line).map((t, tIdx) => {
            if (tIdx < 3) {
                // X, Y, Z coordinates
                t.type = isNumber(t.text) ? 'number' : 'invalid';
            } else if (tIdx < 6) {
                // T/F flags
                t.type = t.text === 'T' || t.text === 'F' ? 'constant' : 'invalid';
            } else {
                t.type = 'comment';
            }
            return t;
        });
    },
    lattVelocitiesStart: (line) => tokenizeConstLine(line, /^[lL]/),
    lattVelocitiesState: (line) => {
        const tokens = tokenizeLine(line);
        tokens.forEach((t, tIdx) => {
            if (tIdx > 0) {
                t.type = 'comment';
            } else if (isInteger(t.text)) {
                t.type = 'number';
            } else {
                t.type = 'invalid';
            }
        });
        return tokens;
    },
    lattVelocitiesVels: tokenizeVector,
    lattVelocitiesLatt: tokenizeVector,
    velocityMode: tokenizeConstLine,
    velocities: tokenizeVector
};

/**
 * Standard tokenizer for vectors (3 numbers).
 */
function tokenizeVector(line: LSPTextLine): Token[] {
    return tokenizeLine(line).map((t, tIdx) => {
        if (tIdx < 3 && isNumber(t.text)) {
            t.type = 'number';
        } else if (tIdx >= 3) {
            t.type = 'comment';
        } else {
            t.type = 'invalid';
        }
        return t;
    });
}

/**
 * Tokenizer for keyword-based lines (e.g. "Direct", "Selective dynamics").
 *
 * @param line - The line to tokenize.
 * @param test - Optional regex to validate the keyword.
 */
function tokenizeConstLine(line: LSPTextLine, test?: RegExp): Token[] {
    const tokens = tokenizeLine(line);
    if (tokens.length > 0) {
        // Check first token against regex
        if (test && !test.test(tokens[0].text)) {
            tokens.forEach((t) => (t.type = 'invalid'));
        } else {
            // Heuristic to distinguish keyword vs comments.
            // VASP comments typically start with #, !, or %.
            let foundComment = false;
            for (const t of tokens) {
                foundComment ||= /^[#!%]/.test(t.text);
                t.type = foundComment ? 'comment' : 'constant';
            }
        }
    }
    return tokens;
}

/**
 * Generic tokenizer that splits a line by whitespace.
 * captures location info (Ranges) for each token.
 */
function tokenizeLine(line: LSPTextLine): Token[] {
    // Regex explanation:
    // ^(\\s*)   -> Match leading whitespace (Group 1)
    // (\\S+)    -> Match a non-whitespace sequence (The token) (Group 2)
    // (.*)$     -> Match everything else (Remainder) (Group 3)
    const matcher = /^(\s*)(\S+)(.*)$/;
    const tokens: Token[] = [];
    let offset = 0;

    // Iteratively peel off tokens from the start of the string
    let matches = line.text.match(matcher);

    while (matches) {
        // matches[1].length is the space before the token
        // matches[2].length is the token itself
        tokens.push({
            range: Range.create(
                line.lineNumber,
                offset + matches[1].length,
                line.lineNumber,
                offset + matches[1].length + matches[2].length
            ),
            text: matches[2]
        });

        // Advance offset past the parsed token
        offset += matches[1].length + matches[2].length;

        // Recurse on the remainder of the line (matches[3])
        matches = matches[3].match(matcher);
    }

    return tokens;
}

/**
 * Helper to sum up atomic counts.
 */
function getNumAtoms(tokens: Token[]): number {
    let numAtoms = 0;
    for (const token of tokens) {
        if (token.type === 'number') {
            numAtoms += +token.text;
        } else {
            break;
        }
    }
    return numAtoms;
}

/**
 * Helper to convert a raw TextDocument into our LSPTextLine abstraction.
 * This effectively pre-processes the document into line objects with range info.
 */
function getDocumentLines(document: TextDocument): LSPTextLine[] {
    const text = document.getText();
    const rawLines = text.split(/\r?\n/);
    return rawLines.map((lineText, i) => {
        const firstNonWhitespace = lineText.search(/\S/);
        return {
            text: lineText,
            lineNumber: i,
            range: Range.create(i, 0, i, lineText.length),
            rangeIncludingLineBreak: Range.create(i, 0, i + 1, 0), // Approximation
            firstNonWhitespaceCharacterIndex: firstNonWhitespace === -1 ? lineText.length : firstNonWhitespace,
            isEmptyOrWhitespace: firstNonWhitespace === -1
        };
    });
}

/**
 * Main Parser Function.
 * Reads the document sequentially and assigns a BlockType to each line based on VASP's fixed structure.
 *
 * @param document - The full text document.
 * @returns Parsed POSCAR document.
 */
export function parsePoscar(document: TextDocument): PoscarDocument {
    const poscarLines: PoscarLine[] = [];
    const lines = getDocumentLines(document);
    const lineCount = lines.length;
    let nextLineIdx = 0;

    /**
     * Consumes one or more lines and attempts to parse them as a specific Type.
     *
     * @param type - Expected block type.
     * @param repeat - Number of lines to consume.
     * @param optionalTest - Additional predicate to validate if a line matches.
     * @returns True if successful, False if validation failed or EOF.
     */
    function processLine(type: PoscarBlockType, repeat?: number, optionalTest?: (tokens: Token[]) => boolean): boolean {
        const myRepeat = repeat ? repeat : 1;
        let isOk = true;
        for (let iter = 0; iter < myRepeat; ++iter) {
            if (lineCount > nextLineIdx) {
                const line = lines[nextLineIdx++];
                const tokens = tokenizers[type](line);

                // If an optional test exists (e.g. checking if species line is strings), fail if it doesn't match
                if (optionalTest && !optionalTest(tokens)) {
                    --nextLineIdx; // Backtrack!
                    isOk = false;
                } else {
                    poscarLines.push({
                        type: type,
                        tokens: tokens,
                        line: line
                    });
                    isOk &&= true;
                }
            }
            // Logic note: why isOk &&= true? It effectively preserves the 'true' state but doesn't handle false correctly?
            // Actually, if isOk becomes false inside, it stays false.
            // If the loop runs and lineCount < nextLineIdx, we don't return false explicitly here?
            isOk &&= true;
        }
        return isOk;
    }

    // --- Parsing State Machine ---

    // 1. Comment Line (Title)
    processLine('comment');

    // 2. Scaling Factor
    processLine('scaling');

    // 3. Lattice Vectors (3x3 Matrix)
    processLine('lattice', 3);

    // 4. Species Names (Optional in VASP 4, Standard in VASP 5)
    // We check if the line contains strings. If not, it might be the atom counts directly.
    processLine('speciesNames', 1, (tokens) => tokens.length > 0 && tokens[0].type === 'string');

    // 5. Number of Atoms
    processLine('numAtoms');

    // Safety check for empty file or partial file
    if (poscarLines.length === 0) return { lines: poscarLines };

    // Calculate total atoms to know how many position lines to read
    const numAtomsLine = poscarLines.find((l) => l.type === 'numAtoms');
    const numAtoms = numAtomsLine ? getNumAtoms(numAtomsLine.tokens) : 0;

    // 6. Selective Dynamics (Optional)
    // Check if line starts with 's' or 'S'.
    const selDyn = processLine('selDynamics', 1, (tokens) => tokens.length > 0 && tokens[0].type === 'constant');

    // 7. Position Mode (Direct/Cartesian)
    processLine('positionMode');

    // 8. Atomic Positions
    // If Selective Dynamics was found, we expect T/F flags.
    processLine(selDyn ? 'positionsSelDyn' : 'positions', numAtoms);

    // 9. Lattice Velocities (Optional, rare)
    if (processLine('lattVelocitiesStart', 1, (tokens) => tokens.length > 0 && tokens[0].type === 'constant')) {
        processLine('lattVelocitiesState');
        processLine('lattVelocitiesVels', 3);
        processLine('lattVelocitiesLatt', 3);
    }

    // 10. Atom Velocities (Optional)
    processLine('velocityMode');
    processLine('velocities', numAtoms);

    return { lines: poscarLines };
}
