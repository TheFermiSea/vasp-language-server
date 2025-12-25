import { FoldingRange, FoldingRangeKind } from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';

export function getFoldingRanges(document: TextDocument): FoldingRange[] {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
    const ranges: FoldingRange[] = [];

    if (document.uri.match(/POSCAR/i) || document.uri.match(/CONTCAR/i)) {
        return getPoscarFoldingRanges(lines);
    } else if (document.uri.match(/INCAR/i)) {
        return getIncarFoldingRanges(lines);
    }

    return ranges;
}

function getPoscarFoldingRanges(lines: string[]): FoldingRange[] {
    const ranges: FoldingRange[] = [];
    if (lines.length < 5) return ranges;

    // 1. Lattice Vectors (Lines 3-5, indices 2-4)
    ranges.push({
        startLine: 2,
        endLine: 4,
        kind: FoldingRangeKind.Region
    });

    // 2. Atomic Coordinates
    // Heuristic: Find where coordinates start
    // Skip Title(0), Scaling(1), Lattice(2,3,4)
    let currentLine = 5;

    // Optional: Species Names (e.g. Fe O)
    if (lines[currentLine] && /^\s*[A-Za-z]/.test(lines[currentLine])) {
        currentLine++;
    }

    // Atom Counts (e.g. 1 1)
    if (lines[currentLine] && /^\s*\d/.test(lines[currentLine])) {
        currentLine++;
    }

    // Optional: Selective Dynamics (starts with S)
    if (lines[currentLine] && /^\s*[Ss]/.test(lines[currentLine])) {
        currentLine++;
    }

    // Coordinate System (Direct/Cartesian)
    if (lines[currentLine] && /^\s*[A-Za-z]/.test(lines[currentLine])) {
        currentLine++;
    }

    // The rest are coordinates
    if (currentLine < lines.length - 1) {
        ranges.push({
            startLine: currentLine,
            endLine: lines.length - 1,
            kind: FoldingRangeKind.Region
        });
    }

    return ranges;
}

function getIncarFoldingRanges(_lines: string[]): FoldingRange[] {
    const ranges: FoldingRange[] = [];
    // future: fold contiguous comment blocks?
    return ranges;
}
