import { Range } from 'vscode-languageserver-types';
import type { BaseToken } from '../types/tokens';

export type { BaseToken } from '../types/tokens';

export function splitLines(text: string): string[] {
    return text.split(/\r?\n/);
}

export function createRange(startLine: number, startChar: number, endLineOrEndChar: number, endChar?: number): Range {
    const endLine = endChar === undefined ? startLine : endLineOrEndChar;
    const finalEndChar = endChar === undefined ? endLineOrEndChar : endChar;
    return Range.create(startLine, startChar, endLine, finalEndChar);
}

export function consumeWhitespace(text: string, offset: number): number {
    const whitespace = text.slice(offset).match(/^\s+/);
    if (whitespace) {
        return offset + whitespace[0].length;
    }
    return offset;
}

export function tokenizeLineByWhitespace(lineText: string, lineNumber: number): BaseToken[] {
    const matcher = /^(\s*)(\S+)(.*)$/;
    const tokens: BaseToken[] = [];
    let offset = 0;

    let matches = lineText.match(matcher);
    while (matches) {
        const leading = matches[1].length;
        const tokenText = matches[2];
        tokens.push({
            range: createRange(lineNumber, offset + leading, offset + leading + tokenText.length),
            text: tokenText
        });
        offset += leading + tokenText.length;
        matches = matches[3].match(matcher);
    }

    return tokens;
}
