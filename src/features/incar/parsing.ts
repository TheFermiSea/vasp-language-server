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
        const continuationIndex = (() => {
            for (let j = line.length - 1; j >= 0; j--) {
                const char = line[j];
                if (char === ' ' || char === '\t') continue;
                return char === '\\' ? j : -1;
            }
            return -1;
        })();
        const lineForParse = continuationIndex >= 0 ? line.slice(0, continuationIndex) : line;
        const isContinuation = continuationIndex >= 0;
        let offset = 0;

        // Check for continuation behavior
        // If line ends with \, it's a continuation.

        while (offset < lineForParse.length) {
            // Skip whitespace
            const whitespace = lineForParse.slice(offset).match(/^\s+/);
            if (whitespace) {
                offset += whitespace[0].length;
                continue;
            }
            if (offset >= lineForParse.length) break;

            const char = lineForParse[offset];
            const remaining = lineForParse.slice(offset);

            // 1. Comment
            if (char === '#' || char === '!') {
                tokens.push({
                    type: 'comment',
                    text: remaining,
                    range: Range.create(i, offset, i, lineForParse.length)
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
                // Inline backslashes are treated as literal characters.
                // A trailing backslash is handled by the line-level continuation check above.
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
                        range: Range.create(i, offset, i, lineForParse.length)
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
    let equalsToken: IncarToken | undefined;
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
                equals: hasEquals ? equalsToken : undefined,
                values: currentValues,
                parsingErrors: diags
            });
        }
        currentTag = null;
        hasEquals = false;
        equalsToken = undefined;
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
                equalsToken = t;
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
