/**
 * CRYSTAL23 .d12 input file parser
 *
 * CRYSTAL23 input files have a multi-section structure:
 * 1. Title line (comment)
 * 2. Geometry section (CRYSTAL/SLAB/POLYMER/MOLECULE/EXTERNAL + structure)
 * 3. Basis set section
 * 4. Hamiltonian/SCF section
 * 5. Properties section (optional)
 *
 * Sections are terminated by END keywords.
 */

import { TextDocument } from 'vscode-languageserver-textdocument';
import { Range, Position } from 'vscode-languageserver-types';
import { createRange, splitLines } from '../../core/parser-utils';
import type { CrystalToken, CrystalTokenType } from '../../types/tokens';
export type { CrystalToken, CrystalTokenType } from '../../types/tokens';

export interface CrystalStatement {
    keyword: CrystalToken;
    values: CrystalToken[];
    range: Range;
}

export type GeometryType = 'CRYSTAL' | 'SLAB' | 'POLYMER' | 'MOLECULE' | 'HELIX' | 'EXTERNAL' | 'UNKNOWN';

export interface GeometrySection {
    type: GeometryType;
    typeToken?: CrystalToken;
    spaceGroup?: number;
    spaceGroupToken?: CrystalToken;
    latticeParams?: number[];
    atomCount?: number;
    atoms: AtomLine[];
    endToken?: CrystalToken;
}

export interface AtomLine {
    atomicNumber: number;
    x: number;
    y: number;
    z: number;
    token: CrystalToken;
}

export interface BasisSetSection {
    entries: BasisSetEntry[];
    endToken?: CrystalToken;
}

export interface BasisSetEntry {
    atomicNumber: number;
    shellCount: number;
    shells: BasisShell[];
}

export interface BasisShell {
    type: string; // S, SP, P, D, F
    primitiveCount: number;
    exponents: number[];
    coefficients: number[];
}

export interface HamiltonianSection {
    statements: CrystalStatement[];
    dftBlock?: DFTBlock;
    endToken?: CrystalToken;
}

export interface DFTBlock {
    functional?: string;
    statements: CrystalStatement[];
}

export interface CrystalDocument {
    title: string;
    titleToken?: CrystalToken;
    geometry: GeometrySection;
    basisSet: BasisSetSection;
    hamiltonian: HamiltonianSection;
    properties?: CrystalStatement[];
    allTokens: CrystalToken[];
    allStatements: CrystalStatement[];
    lines: string[];
    errors: ParseError[];
}

export interface ParseError {
    message: string;
    range: Range;
    severity: 'error' | 'warning';
}

// Geometry type keywords
const GEOMETRY_TYPES = ['CRYSTAL', 'SLAB', 'POLYMER', 'MOLECULE', 'HELIX', 'EXTERNAL'];

// Keywords that start blocks
const BLOCK_START_KEYWORDS = ['DFT', 'OPTGEOM', 'FREQCALC', 'ELASTCON', 'EOS', 'BAND', 'DOSS', 'BASISSET'];

// Keywords that are standalone
const STANDALONE_KEYWORDS = [
    // SCF
    'SHRINK',
    'TOLINTEG',
    'FMIXING',
    'BROYDEN',
    'DIIS',
    'ANDERSON',
    'MAXCYCLE',
    'LEVSHIFT',
    'SPINLOCK',
    'TOLDEE',
    'TOLDEP',
    'TOLDEG',
    // DFT functionals
    'B3LYP',
    'PBE',
    'PBESOL',
    'HSE06',
    'M06',
    'PW91',
    'LDA',
    'LSDA',
    'SPIN',
    'HFEXCHANGE',
    'XLGRID',
    'XXLGRID',
    // Geometry optimization
    'FULLOPTG',
    'ATOMONLY',
    'CELLONLY',
    'CVOLOPT',
    'INTRANS',
    // Properties
    'ECHG',
    'POTC',
    // Misc
    'ENDSCF',
    'ENDOPT',
    'ENDDFT',
    'NOCOMM'
];

/**
 * Parse a CRYSTAL23 .d12 input file into sections and tokens.
 *
 * @param document - LSP text document for a CRYSTAL23 input file.
 * @returns Parsed CRYSTAL document with tokens, statements, and errors.
 */
export function parseCrystal(document: TextDocument): CrystalDocument {
    const text = document.getText();
    const lines = splitLines(text);
    const allTokens: CrystalToken[] = [];
    const allStatements: CrystalStatement[] = [];
    const errors: ParseError[] = [];

    // Minimum file validation
    if (lines.length < 2) {
        errors.push({
            message: `CRYSTAL23 input file is too short (${lines.length} line${lines.length === 1 ? '' : 's'}). A valid .d12 file requires at minimum:\n  Line 1: Title/comment\n  Line 2: Geometry type (CRYSTAL, SLAB, POLYMER, MOLECULE)\n  Plus: Space group, lattice parameters, atom coordinates, basis set, and SCF settings.`,
            range: createRange(0, 0, 0, lines[0]?.length || 0),
            severity: 'error'
        });
    }

    // Parse title (first line)
    const titleToken = createToken('title', lines[0] || '', 0, 0, lines[0]?.length || 0);
    allTokens.push(titleToken);

    // Warn if title is empty
    if (!lines[0] || lines[0].trim().length === 0) {
        errors.push({
            message:
                'Line 1: Empty title line. While valid, it is recommended to include a descriptive title for your calculation.',
            range: createRange(0, 0, 0, 0),
            severity: 'warning'
        });
    }

    // Initialize sections
    const geometry: GeometrySection = {
        type: 'UNKNOWN',
        atoms: []
    };
    const basisSet: BasisSetSection = { entries: [] };
    const hamiltonian: HamiltonianSection = { statements: [] };

    let currentLine = 1;

    // Parse geometry type (line 2)
    if (lines.length > 1) {
        const geomLine = lines[1].trim().toUpperCase();
        for (const gt of GEOMETRY_TYPES) {
            if (geomLine.startsWith(gt)) {
                geometry.type = gt as GeometryType;
                geometry.typeToken = createToken('geometry-type', gt, 1, 0, gt.length);
                allTokens.push(geometry.typeToken);
                break;
            }
        }
        if (geometry.type === 'UNKNOWN') {
            errors.push({
                message: `Line 2: Unknown geometry type '${geomLine}'. Expected one of:\n  - CRYSTAL: 3D periodic structure\n  - SLAB: 2D periodic (surface)\n  - POLYMER: 1D periodic (chain)\n  - MOLECULE: Non-periodic (cluster)\n  - HELIX: Helical symmetry\n  - EXTERNAL: Read geometry from external file`,
                range: createRange(1, 0, 1, lines[1].length),
                severity: 'error'
            });
        }
        currentLine = 2;
    } else {
        errors.push({
            message:
                'Line 2: Missing geometry type declaration. A CRYSTAL23 input file must specify the system dimensionality (CRYSTAL, SLAB, POLYMER, MOLECULE, etc.).',
            range: createRange(1, 0, 1, 0),
            severity: 'error'
        });
    }

    // Parse remaining lines based on context
    let inBasisSet = false;
    let inDFTBlock = false;

    for (let i = currentLine; i < lines.length; i++) {
        const line = lines[i];
        const trimmed = line.trim();
        const upper = trimmed.toUpperCase();

        // Skip empty lines
        if (!trimmed) continue;

        // Check for END keyword
        if (upper === 'END' || upper === 'ENDG' || upper === 'ENDB' || upper === 'ENDSCF' || upper === 'ENDDFT') {
            const endToken = createToken('end', trimmed, i, 0, trimmed.length);
            allTokens.push(endToken);

            if (inDFTBlock) {
                inDFTBlock = false;
                hamiltonian.dftBlock!.statements = hamiltonian.statements.slice();
            }
            if (inBasisSet) {
                inBasisSet = false;
                basisSet.endToken = endToken;
            }
            continue;
        }

        // Check for comment lines
        if (trimmed.startsWith('#')) {
            const commentToken = createToken('comment', trimmed, i, 0, trimmed.length);
            allTokens.push(commentToken);
            continue;
        }

        // Check for block-starting keywords
        if (BLOCK_START_KEYWORDS.includes(upper.split(/\s+/)[0])) {
            const keyword = upper.split(/\s+/)[0];
            const keywordToken = createToken('block-start', keyword, i, 0, keyword.length);
            allTokens.push(keywordToken);

            if (keyword === 'DFT') {
                inDFTBlock = true;
                hamiltonian.dftBlock = { statements: [] };
            } else if (keyword === 'BASISSET') {
                inBasisSet = true;
            }

            const statement: CrystalStatement = {
                keyword: keywordToken,
                values: [],
                range: createRange(i, 0, i, line.length)
            };
            allStatements.push(statement);
            hamiltonian.statements.push(statement);
            continue;
        }

        // Check for standalone keywords
        const words = upper.split(/\s+/);
        if (STANDALONE_KEYWORDS.includes(words[0])) {
            const keyword = words[0];
            const keywordToken = createToken('keyword', keyword, i, 0, keyword.length);
            allTokens.push(keywordToken);

            // Parse values after keyword
            const values: CrystalToken[] = [];
            let pos = keyword.length;
            for (let j = 1; j < words.length; j++) {
                const startPos = line.indexOf(words[j], pos);
                const valueToken = createToken('number', words[j], i, startPos, startPos + words[j].length);
                values.push(valueToken);
                allTokens.push(valueToken);
                pos = startPos + words[j].length;
            }

            const statement: CrystalStatement = {
                keyword: keywordToken,
                values,
                range: createRange(i, 0, i, line.length)
            };
            allStatements.push(statement);
            hamiltonian.statements.push(statement);
            continue;
        }

        // Check if line looks like a number sequence (lattice params, atom coords, etc.)
        if (/^[\d.\s\-+eE]+$/.test(trimmed)) {
            const numToken = createToken('number', trimmed, i, 0, trimmed.length);
            allTokens.push(numToken);

            // Validate numeric values
            const numParts = trimmed.split(/\s+/).filter((p) => p.length > 0);
            for (let j = 0; j < numParts.length; j++) {
                const val = Number(numParts[j]);
                if (isNaN(val)) {
                    errors.push({
                        message: `Line ${i + 1}: Invalid numeric value '${numParts[j]}' at position ${j + 1}. Expected a valid number (integer or decimal, with optional scientific notation like 1.23e-4).`,
                        range: createRange(i, 0, i, trimmed.length),
                        severity: 'error'
                    });
                }
            }
            continue;
        }

        // Unknown token - provide context-aware suggestions
        const unknownToken = createToken('unknown', trimmed, i, 0, trimmed.length);
        allTokens.push(unknownToken);

        // Try to provide helpful suggestions for common mistakes
        const upperTrimmed = trimmed.toUpperCase();
        let suggestion = '';

        // Check for common typos or similar keywords
        const similarKeywords = [...BLOCK_START_KEYWORDS, ...STANDALONE_KEYWORDS].filter(
            (kw) => kw.includes(upperTrimmed.slice(0, 3)) || upperTrimmed.includes(kw.slice(0, 3))
        );
        if (similarKeywords.length > 0) {
            suggestion = ` Did you mean: ${similarKeywords.slice(0, 3).join(', ')}?`;
        }

        errors.push({
            message: `Line ${i + 1}: Unrecognized keyword or data '${trimmed}'.${suggestion} If this is numeric data, ensure it contains only numbers, spaces, and decimal points.`,
            range: createRange(i, 0, i, trimmed.length),
            severity: 'warning'
        });
    }

    return {
        title: lines[0] || '',
        titleToken,
        geometry,
        basisSet,
        hamiltonian,
        allTokens,
        allStatements,
        lines,
        errors
    };
}

function createToken(
    type: CrystalTokenType,
    text: string,
    line: number,
    startChar: number,
    endChar: number
): CrystalToken {
    return {
        type,
        text,
        range: createRange(line, startChar, line, endChar),
        lineNumber: line
    };
}

/**
 * Get the token at a given position
 */
/**
 * Locate the token at a given position in the parsed CRYSTAL document.
 *
 * @param document - Parsed CRYSTAL document.
 * @param position - LSP position to look up.
 * @returns The token at the position, or undefined if none matches.
 */
export function getTokenAtPosition(document: CrystalDocument, position: Position): CrystalToken | undefined {
    return document.allTokens.find((token) => isPositionInRange(position, token.range));
}

function isPositionInRange(position: Position, range: Range): boolean {
    if (position.line < range.start.line || position.line > range.end.line) return false;
    if (position.line === range.start.line && position.character < range.start.character) return false;
    if (position.line === range.end.line && position.character > range.end.character) return false;
    return true;
}
