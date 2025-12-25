import { TextDocument } from 'vscode-languageserver-textdocument';
import { Range, Diagnostic, DiagnosticSeverity } from 'vscode-languageserver-types';

/**
 * Types of Tokens in INCAR.
 */
export type IncarTokenType = 'tag' | 'equals' | 'value' | 'comment' | 'semicolon' | 'eol' | 'invalid';

export interface IncarToken {
    type: IncarTokenType;
    text: string;
    range: Range;
}

/**
 * Represents a parsed "Tag = Value" statement.
 */
export interface IncarStatement {
    tag: IncarToken;
    /** The equals sign token (optional/missing if syntax error) */
    equals?: IncarToken;
    /** All tokens that make up the value (could be multiple if array or complex string) */
    values: IncarToken[];
    /** Any validation errors specific to parsing (e.g., missing equals) */
    parsingErrors: Diagnostic[];
}

/**
 * Represents the full parsed result of an INCAR file.
 */
export interface IncarDocument {
    statements: IncarStatement[];
    allTokens: IncarToken[];
}

/**
 * Parses an INCAR document.
 */
export function parseIncar(document: TextDocument): IncarDocument {
    const text = document.getText();
    const tokens = tokenizeIncar(text, document);
    const statements = groupTokensIntoStatements(tokens);
    return { statements, allTokens: tokens };
}

/**
 * Tokenizer for INCAR.
 */
function tokenizeIncar(text: string, _document: TextDocument): IncarToken[] {
    const tokens: IncarToken[] = [];
    const lines = text.split(/\r?\n/);

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i];
        let offset = 0;
        let isContinuation = false;

        // Check for continuation behavior
        // If line ends with \, it's a continuation.

        while (offset < line.length) {
            // Skip whitespace
            const whitespace = line.slice(offset).match(/^\s+/);
            if (whitespace) {
                offset += whitespace[0].length;
                continue;
            }
            if (offset >= line.length) break;

            const char = line[offset];
            const remaining = line.slice(offset);

            // 1. Comment
            if (char === '#' || char === '!') {
                tokens.push({
                    type: 'comment',
                    text: remaining,
                    range: Range.create(i, offset, i, line.length)
                });
                break; // Rest of line is comment
            }

            // 2. Semicolon
            if (char === ';') {
                tokens.push({
                    type: 'semicolon',
                    text: ';',
                    range: Range.create(i, offset, i, offset + 1)
                });
                offset++;
                continue;
            }

            // 3. Equals
            if (char === '=') {
                tokens.push({
                    type: 'equals',
                    text: '=',
                    range: Range.create(i, offset, i, offset + 1)
                });
                offset++;
                continue;
            }

            // 4. Line Continuation (\)
            if (char === '\\') {
                // If followed by optional whitespace and then newline?
                // We just see if it's the last significant char.
                // For now, treat as token, parser state logic handles "no EOL".
                isContinuation = true;
                offset++;
                continue;
            }

            // 5. String
            if (char === '"' || char === "'") {
                const quoteRegex = char === '"' ? /^"[^"]*"/ : /^'[^']*'/;
                const match = remaining.match(quoteRegex);
                if (match) {
                    tokens.push({
                        type: 'value',
                        text: match[0],
                        range: Range.create(i, offset, i, offset + match[0].length)
                    });
                    offset += match[0].length;
                    continue;
                } else {
                    // Unclosed string
                    tokens.push({
                        type: 'invalid',
                        text: remaining,
                        range: Range.create(i, offset, i, line.length)
                    });
                    break;
                }
            }

            // 6. Word
            const wordMatch = remaining.match(/^[^=\s;#!\\]+/);
            if (wordMatch) {
                const word = wordMatch[0];
                tokens.push({
                    type: 'value',
                    text: word,
                    range: Range.create(i, offset, i, offset + word.length)
                });
                offset += word.length;
            } else {
                offset++;
            }
        }

        if (!isContinuation) {
            tokens.push({
                type: 'eol',
                text: '\n',
                range: Range.create(i, line.length, i, line.length)
            });
        }
    }
    return tokens;
}

/**
 * Contextual grouping: Tag = Value
 */
function groupTokensIntoStatements(tokens: IncarToken[]): IncarStatement[] {
    const statements: IncarStatement[] = [];
    let currentTag: IncarToken | null = null;
    let hasEquals = false;
    let currentValues: IncarToken[] = [];

    // Helper to flush current statement
    const flush = () => {
        if (currentTag) {
            const diags: Diagnostic[] = [];
            if (!hasEquals) {
                diags.push({
                    message: `Expected '=' after tag '${currentTag.text}'.`,
                    range: currentTag.range,
                    severity: DiagnosticSeverity.Error
                });
            } else if (currentValues.length === 0) {
                diags.push({
                    message: `Tag '${currentTag.text}' has no value assigned.`,
                    range: currentTag.range,
                    severity: DiagnosticSeverity.Warning
                });
            }

            statements.push({
                tag: { ...currentTag, type: 'tag' },
                equals: hasEquals ? { type: 'equals', text: '=', range: Range.create(0, 0, 0, 0) } : undefined, // Placeholder range need fix? No, just existence check mostly.
                values: currentValues,
                parsingErrors: diags
            });
        }
        currentTag = null;
        hasEquals = false;
        currentValues = [];
    };

    let i = 0;
    while (i < tokens.length) {
        const t = tokens[i];

        if (t.type === 'comment') {
            // Comments are ignored
            i++;
            continue;
        }

        if (t.type === 'semicolon' || t.type === 'eol') {
            flush();
            i++;
            continue;
        }

        if (t.type === 'equals') {
            if (currentTag) {
                hasEquals = true;
            }
            i++;
            continue;
        }

        if (t.type === 'value' || t.type === 'invalid') {
            if (!currentTag) {
                currentTag = t;
            } else {
                if (!hasEquals) {
                    // Missing equals logic?
                    // Try to be smart? No, strictly require equals.
                    // If we see TAG KEYWORD without equals, it is messy.
                    // But assume current is value and let linter complain about missing equals.
                    currentValues.push(t);
                } else {
                    currentValues.push(t);
                }
            }
            i++;
        } else {
            i++;
        }
    }
    flush();

    return statements;
}
