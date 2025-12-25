import { TextDocument } from 'vscode-languageserver-textdocument';
import { TextEdit, Range, Position } from 'vscode-languageserver-types';

export function formatIncar(document: TextDocument): TextEdit[] {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
    const edits: TextEdit[] = [];

    // 1. Analyze lines to find max key length
    let maxKeyLength = 0;
    const splitLines: { key: string; value: string; originalIndex: number; isAssignment: boolean }[] = [];

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i];
        const trimmed = line.trim();

        // Skip comments and empty lines for alignment calculation
        if (trimmed.length === 0 || trimmed.startsWith('#') || trimmed.startsWith('!')) {
            splitLines.push({ key: line, value: '', originalIndex: i, isAssignment: false });
            continue;
        }

        const eqIndex = line.indexOf('='); // Naive check, assumes first = is operator
        if (eqIndex !== -1) {
            const key = line.substring(0, eqIndex).trim();
            // Value includes comment if present
            const value = line.substring(eqIndex + 1).trim();

            if (key.length > maxKeyLength) {
                maxKeyLength = key.length;
            }
            splitLines.push({ key, value, originalIndex: i, isAssignment: true });
        } else {
            splitLines.push({ key: line, value: '', originalIndex: i, isAssignment: false });
        }
    }

    // 2. Generate Edits
    // We replace the entire document content to be safe and simple,
    // or we can generate per-line edits. Per-line is better for minimal diffs.
    // However, aligning affects many lines.
    // Let's replace line by line if different.

    for (const item of splitLines) {
        if (item.isAssignment) {
            // Pad the key
            // Output format: KEY   = VALUE
            const paddedKey = item.key.padEnd(maxKeyLength, ' ');
            const newText = `${paddedKey} = ${item.value}`;

            if (newText !== lines[item.originalIndex]) {
                const range = Range.create(
                    Position.create(item.originalIndex, 0),
                    Position.create(item.originalIndex, lines[item.originalIndex].length)
                );
                edits.push(TextEdit.replace(range, newText));
            }
        }
    }

    return edits;
}
