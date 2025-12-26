import { DocumentSymbol, SymbolKind, Range } from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { IncarDocument, parseIncar } from './parsing';

/**
 * Generates Document Symbols for INCAR files.
 * Maps every TAG = VALUE pair to a Variable symbol.
 */
export function getIncarSymbols(document: TextDocument, parsed?: IncarDocument): DocumentSymbol[] {
    const symbols: DocumentSymbol[] = [];
    const myParsed = parsed || parseIncar(document);

    for (const stmt of myParsed.statements) {
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
