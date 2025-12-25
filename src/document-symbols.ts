import { DocumentSymbol, SymbolKind, Range } from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { parseIncar } from './incar-parsing';

/**
 * Generates Document Symbols for INCAR files.
 * Maps every TAG = VALUE pair to a Variable symbol.
 */
export function getIncarSymbols(document: TextDocument): DocumentSymbol[] {
    const symbols: DocumentSymbol[] = [];
    const parsed = parseIncar(document);

    for (const stmt of parsed.statements) {
        // Tag name as the symbol name
        const name = stmt.tag.text;

        // Use the whole range of the statement (Tag ... Value) for the symbol
        const range = Range.create(
            stmt.tag.range.start,
            stmt.values.length > 0 ? stmt.values[stmt.values.length - 1].range.end : stmt.tag.range.end
        );

        // Selection range is just the Tag itself
        const selectionRange = stmt.tag.range;

        // Detail shows the value
        const detail = stmt.values.map((v) => v.text).join(' ');

        symbols.push(DocumentSymbol.create(name, detail, SymbolKind.Variable, range, selectionRange));
    }

    return symbols;
}

/**
 * Generates Document Symbols for POSCAR files.
 * Uses heuristics to identify sections (Lattice, Species, Coordinates).
 */
export function getPoscarSymbols(document: TextDocument): DocumentSymbol[] {
    const symbols: DocumentSymbol[] = [];
    const text = document.getText();
    const lines = text.split(/\r?\n/);

    if (lines.length < 5) return symbols;

    // 1. Comment / Title (Line 1)
    symbols.push(
        DocumentSymbol.create(
            'Title/Comment',
            lines[0].trim(),
            SymbolKind.String,
            Range.create(0, 0, 0, lines[0].length),
            Range.create(0, 0, 0, lines[0].length)
        )
    );

    // 2. Scaling Factor (Line 2)
    symbols.push(
        DocumentSymbol.create(
            'Scaling Factor',
            lines[1].trim(),
            SymbolKind.Number,
            Range.create(1, 0, 1, lines[1].length),
            Range.create(1, 0, 1, lines[1].length)
        )
    );

    // 3. Lattice Vectors (Lines 3-5)
    symbols.push(
        DocumentSymbol.create(
            'Lattice Vectors',
            'a1, a2, a3',
            SymbolKind.Array,
            Range.create(2, 0, 4, lines[4].length),
            Range.create(2, 0, 2, lines[2].length)
        )
    );

    // Heuristics for Species and Atoms
    // Line 6 might be Species names (Letters) or Count (Numbers)
    let currentLine = 5;
    if (lines[currentLine] && /^\s*[A-Za-z]/.test(lines[currentLine])) {
        symbols.push(
            DocumentSymbol.create(
                'Species Names',
                lines[currentLine].trim(),
                SymbolKind.Field,
                Range.create(currentLine, 0, currentLine, lines[currentLine].length),
                Range.create(currentLine, 0, currentLine, lines[currentLine].length)
            )
        );
        currentLine++;
    }

    if (lines[currentLine] && /^\s*\d/.test(lines[currentLine])) {
        symbols.push(
            DocumentSymbol.create(
                'Atom Counts',
                lines[currentLine].trim(),
                SymbolKind.Number,
                Range.create(currentLine, 0, currentLine, lines[currentLine].length),
                Range.create(currentLine, 0, currentLine, lines[currentLine].length)
            )
        );
        currentLine++;
    }

    // Direct/Cartesian
    if (lines[currentLine]) {
        symbols.push(
            DocumentSymbol.create(
                'Coordinate System',
                lines[currentLine].trim(),
                SymbolKind.EnumMember,
                Range.create(currentLine, 0, currentLine, lines[currentLine].length),
                Range.create(currentLine, 0, currentLine, lines[currentLine].length)
            )
        );
        currentLine++;
    }

    // Coordinates (remaining lines)
    if (currentLine < lines.length) {
        symbols.push(
            DocumentSymbol.create(
                'Atomic Coordinates',
                `${lines.length - currentLine} entries`,
                SymbolKind.Struct,
                Range.create(currentLine, 0, lines.length - 1, lines[lines.length - 1].length),
                Range.create(currentLine, 0, currentLine, lines[currentLine].length)
            )
        );
    }

    return symbols;
}
