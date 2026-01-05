import { FoldingRange, FoldingRangeKind } from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';

/**
 * Compute folding ranges for INCAR section headers.
 *
 * @param document - LSP text document for an INCAR file.
 * @returns Folding ranges representing comment-delimited sections.
 */
export function getIncarFoldingRanges(document: TextDocument): FoldingRange[] {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
    const ranges: FoldingRange[] = [];

    // Region-based folding on comments like "# --- Section ---"
    // Also basic indentation folding or block folding.
    // Since INCAR is flat, we'll try to fold "sections" defined by header comments.
    // Strategy: A header comment starts a section. The section ends at the next header comment or EOF.

    const headerRegex = /^#\s*([A-Z].*|-+.*)/; // Matches "# ---" or "# Section"
    let currentStartLine: number | null = null;

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();

        if (headerRegex.test(line)) {
            // Close previous section if open
            if (currentStartLine !== null && i > currentStartLine + 1) {
                ranges.push(
                    FoldingRange.create(currentStartLine, i - 1, undefined, undefined, FoldingRangeKind.Region)
                );
            }
            currentStartLine = i;
        }
    }

    // Close last section
    if (currentStartLine !== null && lines.length - 1 > currentStartLine) {
        ranges.push(
            FoldingRange.create(currentStartLine, lines.length - 1, undefined, undefined, FoldingRangeKind.Region)
        );
    }

    return ranges;
}
