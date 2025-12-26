import { SemanticTokensBuilder } from 'vscode-languageserver/node';
import { tokenTypeIndices } from '../semantic-tokens-legend';

export function getKpointsSemanticTokens(text: string): { data: number[] } {
    const builder = new SemanticTokensBuilder();
    const lines = text.split(/\r?\n/);

    lines.forEach((lineText, lineIdx) => {
        // Simple regex-based tokenization since KPOINTS is line-based
        const words = lineText.match(/(\S+)/g);
        if (!words) return;

        let currentOffset = 0;

        // Find positions of words
        // Note: regex match doesn't give indices easily in loop, so we re-scan
        const wordMatches = lineText.matchAll(/(\S+)/g);

        for (const match of wordMatches) {
            const word = match[0];
            const startChar = match.index || 0;

            let typeIdx: number | undefined;

            if (lineIdx === 0) {
                typeIdx = tokenTypeIndices.get('comment');
            } else if (lineIdx === 2) {
                // Mode line
                typeIdx = tokenTypeIndices.get('keyword');
            } else {
                if (!isNaN(+word)) {
                    typeIdx = tokenTypeIndices.get('number');
                } else {
                    typeIdx = tokenTypeIndices.get('string'); // Fallback for comments on data lines?
                }
            }

            // Override: Comments starting with !, #
            if (/^[!#]/.test(word)) {
                typeIdx = tokenTypeIndices.get('comment');
                // Color rest of logic as comment? Simplified for now.
            }

            if (typeIdx !== undefined) {
                builder.push(lineIdx, startChar, word.length, typeIdx, 0);
            }
        }
    });

    return builder.build();
}
