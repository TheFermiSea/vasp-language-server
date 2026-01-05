import { SemanticTokensBuilder } from 'vscode-languageserver/node';
import { tokenTypeIndices } from '../semantic-tokens-legend';

/**
 * Build semantic tokens for a KPOINTS document.
 *
 * @param text - Full document text.
 * @returns Semantic tokens payload for the LSP response.
 */
export function getKpointsSemanticTokens(text: string): { data: number[] } {
    const builder = new SemanticTokensBuilder();
    const lines = text.split(/\r?\n/);

    const numberType = tokenTypeIndices.get('number');
    const keywordType = tokenTypeIndices.get('keyword');
    const commentType = tokenTypeIndices.get('comment');
    const stringType = tokenTypeIndices.get('string');

    const addLineTokens = (lineText: string | undefined, lineIdx: number, mode: 'auto' | 'number' | 'keyword') => {
        if (!lineText) return;
        const commentStart = lineText.search(/[!#]/);
        const hasComment = commentStart >= 0;
        const dataPart = hasComment ? lineText.slice(0, commentStart) : lineText;
        const commentPart = hasComment ? lineText.slice(commentStart) : '';

        for (const match of dataPart.matchAll(/\S+/g)) {
            const word = match[0];
            const startChar = match.index ?? 0;
            let typeIdx: number | undefined;

            if (mode === 'keyword') {
                typeIdx = keywordType;
            } else if (mode === 'number') {
                typeIdx = numberType ?? stringType;
            } else {
                typeIdx = !Number.isNaN(Number(word)) ? numberType : stringType;
            }

            if (typeIdx !== undefined) {
                builder.push(lineIdx, startChar, word.length, typeIdx, 0);
            }
        }

        if (commentPart.trim().length > 0 && commentType !== undefined) {
            builder.push(lineIdx, commentStart, commentPart.length, commentType, 0);
        }
    };

    // Line 1: comment/description
    if (lines[0] && commentType !== undefined) {
        builder.push(0, 0, lines[0].length, commentType, 0);
    }

    // Line 2: number of k-points
    addLineTokens(lines[1], 1, 'number');

    // Line 3: mode
    addLineTokens(lines[2], 2, 'keyword');

    const rawCount = lines[1]?.trim().split(/\s+/)[0] ?? '';
    const numKpoints = Number.parseInt(rawCount, 10);
    const modeText = (lines[2] ?? '').trim().toLowerCase();
    const isExplicit = Number.isFinite(numKpoints) && numKpoints > 0;
    const isLineMode = modeText.startsWith('l');

    if (isExplicit) {
        // Explicit k-points: line 4 is coordinate system or line-mode qualifier.
        addLineTokens(lines[3], 3, 'keyword');
        const dataStart = 4;
        for (let i = dataStart; i < lines.length; i++) {
            addLineTokens(lines[i], i, 'auto');
        }
        return builder.build();
    }

    // Automatic modes: line 4 is grid or length, line 5 optional shift
    addLineTokens(lines[3], 3, 'auto');
    addLineTokens(lines[4], 4, 'auto');

    if (isLineMode) {
        // Line-mode can still include labeled point lists after the header.
        for (let i = 5; i < lines.length; i++) {
            addLineTokens(lines[i], i, 'auto');
        }
    }

    return builder.build();
}
