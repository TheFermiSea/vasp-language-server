import { FoldingRange, FoldingRangeKind } from 'vscode-languageserver/node';
import { PoscarLine } from './parsing';

export function getFoldingRanges(lines: string[], poscarData?: PoscarLine[]): FoldingRange[] {
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
    // If we have poscarData, we can be more precise.
    let currentLine = 5;

    if (poscarData) {
        const firstPos = poscarData.find((l) => l.type === 'positions' || l.type === 'positionsSelDyn');
        if (firstPos) {
            currentLine = firstPos.line.lineNumber;
        }
    } else {
        // Skip Title(0), Scaling(1), Lattice(2,3,4)

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
