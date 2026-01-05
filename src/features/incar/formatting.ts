import { TextDocument } from 'vscode-languageserver-textdocument';
import { TextEdit, Range, Position } from 'vscode-languageserver-types';
import { parseIncar } from './parsing';

/**
 * Represents a parsed line for formatting purposes.
 */
interface FormattableLine {
    /** The tag/key name */
    key: string;
    /** The value(s) as a string (may include trailing comment) */
    value: string;
    /** Original line number (0-indexed) */
    lineNumber: number;
    /** Whether this is a valid assignment (tag = value) */
    isAssignment: boolean;
    /** End column of the last token on this line (for range calculation) */
    lineEndColumn: number;
}

/**
 * Formats an INCAR document by aligning all tag = value assignments.
 * Uses the tokenized AST from parseIncar for robust parsing that correctly
 * handles edge cases like comments containing '=' characters.
 *
 * @param document - LSP text document for an INCAR file.
 * @returns Text edits that align assignments without changing content.
 */
export function formatIncar(document: TextDocument): TextEdit[] {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
    const edits: TextEdit[] = [];

    // Parse document using the robust tokenizer
    const parsed = parseIncar(document);

    // Build a map of line number -> formattable line info from statements
    const lineMap = new Map<number, FormattableLine>();

    for (const stmt of parsed.statements) {
        // Get line number from tag token
        const lineNumber = stmt.tag.range.start.line;

        // Build value string from value tokens
        const valueStr = stmt.values.map((v) => v.text).join(' ');

        // Find trailing comment on this line (if any)
        const trailingComment = findTrailingComment(parsed.allTokens, lineNumber);
        const fullValue = trailingComment ? `${valueStr} ${trailingComment}`.trim() : valueStr;

        // Calculate line end for range
        const lineEndColumn = lines[lineNumber]?.length ?? 0;

        lineMap.set(lineNumber, {
            key: stmt.tag.text,
            value: fullValue,
            lineNumber,
            isAssignment: stmt.equals !== undefined,
            lineEndColumn
        });
    }

    // Calculate max key length from all valid assignments
    let maxKeyLength = 0;
    for (const item of lineMap.values()) {
        if (item.isAssignment && item.key.length > maxKeyLength) {
            maxKeyLength = item.key.length;
        }
    }

    // Generate edits for each assignment line
    for (const item of lineMap.values()) {
        if (item.isAssignment) {
            const paddedKey = item.key.padEnd(maxKeyLength, ' ');
            const newText = `${paddedKey} = ${item.value}`;
            const originalLine = lines[item.lineNumber];

            if (newText !== originalLine) {
                const range = Range.create(
                    Position.create(item.lineNumber, 0),
                    Position.create(item.lineNumber, originalLine.length)
                );
                edits.push(TextEdit.replace(range, newText));
            }
        }
    }

    return edits;
}

/**
 * Finds a trailing comment token on the specified line.
 * Returns the comment text if found, undefined otherwise.
 */
function findTrailingComment(
    tokens: { type: string; text: string; range: Range }[],
    lineNumber: number
): string | undefined {
    for (const token of tokens) {
        if (token.type === 'comment' && token.range.start.line === lineNumber) {
            return token.text;
        }
    }
    return undefined;
}
