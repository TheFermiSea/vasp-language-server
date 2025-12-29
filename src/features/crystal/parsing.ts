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

export type CrystalTokenType =
    | 'title'
    | 'geometry-type'
    | 'space-group'
    | 'lattice-param'
    | 'atom-count'
    | 'atom-line'
    | 'keyword'
    | 'block-start'
    | 'block-end'
    | 'basis-header'
    | 'basis-data'
    | 'number'
    | 'comment'
    | 'end'
    | 'unknown';

export interface CrystalToken {
    type: CrystalTokenType;
    text: string;
    range: Range;
    lineNumber: number;
}

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
 * Parse a CRYSTAL23 .d12 input file
 */
export function parseCrystal(document: TextDocument): CrystalDocument {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
    const allTokens: CrystalToken[] = [];
    const allStatements: CrystalStatement[] = [];
    const errors: ParseError[] = [];

    // Parse title (first line)
    const titleToken = createToken('title', lines[0] || '', 0, 0, lines[0]?.length || 0);
    allTokens.push(titleToken);

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
                message: `Unknown geometry type: ${geomLine}. Expected: ${GEOMETRY_TYPES.join(', ')}`,
                range: createRange(1, 0, 1, lines[1].length),
                severity: 'error'
            });
        }
        currentLine = 2;
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
        if (/^[\d.\s-]+$/.test(trimmed)) {
            const numToken = createToken('number', trimmed, i, 0, trimmed.length);
            allTokens.push(numToken);
            continue;
        }

        // Unknown token
        const unknownToken = createToken('unknown', trimmed, i, 0, trimmed.length);
        allTokens.push(unknownToken);
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

function createRange(startLine: number, startChar: number, endLine: number, endChar: number): Range {
    return {
        start: Position.create(startLine, startChar),
        end: Position.create(endLine, endChar)
    };
}

/**
 * Get the token at a given position
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
