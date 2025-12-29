import { CompletionItem, CompletionItemKind, Position } from 'vscode-languageserver/node';

export function getKpointsCompletions(text: string, position: Position): CompletionItem[] {
    // Only provide completions for the 3rd line (Index 2), which is the Mode.
    // Line 0: Comment
    // Line 1: Number of kpoints (0 for auto)
    // Line 2: Mode
    if (position.line !== 2) return [];

    return [
        {
            label: 'Automatic',
            kind: CompletionItemKind.Keyword,
            detail: 'KPOINTS Generation Mode',
            documentation: 'Fully automatic mesh generation. Requires 0 number of k-points.'
        },
        {
            label: 'Gamma',
            kind: CompletionItemKind.Keyword,
            detail: 'KPOINTS Generation Mode',
            documentation: 'Gamma-centered grid. Recommended for Hexagonal lattices.'
        },
        {
            label: 'Monkhorst-Pack',
            kind: CompletionItemKind.Keyword,
            detail: 'KPOINTS Generation Mode',
            documentation: 'Standard mesh generation shifting away from Gamma.'
        },
        {
            label: 'Line-Mode',
            kind: CompletionItemKind.Keyword,
            detail: 'KPOINTS Generation Mode',
            documentation: 'For Band Structure calculations along high-symmetry lines.'
        },
        // Snippets
        {
            label: 'Snippet: Monkhorst-Pack 4x4x4',
            kind: CompletionItemKind.Snippet,
            detail: 'MP Grid Template',
            insertText: 'Monkhorst-Pack\n${1:4} ${2:4} ${3:4}\n0 0 0',
            insertTextFormat: 2 // Snippet
        },
        {
            label: 'Snippet: Gamma 4x4x4',
            kind: CompletionItemKind.Snippet,
            detail: 'Gamma Grid Template',
            insertText: 'Gamma\n${1:4} ${2:4} ${3:4}\n0 0 0',
            insertTextFormat: 2 // Snippet
        }
    ];
}
