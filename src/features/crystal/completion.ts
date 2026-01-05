/**
 * CRYSTAL23 keyword completion
 */

import { CompletionItem, CompletionItemKind, InsertTextFormat } from 'vscode-languageserver-types';
import { CRYSTAL_TAGS, CrystalTagDefinition } from '../../data/crystal-tags';

/**
 * Get completion items for CRYSTAL23 keywords.
 *
 * @returns Completion items derived from the CRYSTAL tag database.
 */
export function getCrystalCompletions(): CompletionItem[] {
    const items: CompletionItem[] = [];

    for (const [keyword, def] of Object.entries(CRYSTAL_TAGS)) {
        items.push(createCompletionItem(keyword, def));
    }

    return items;
}

function createCompletionItem(keyword: string, def: CrystalTagDefinition): CompletionItem {
    const item: CompletionItem = {
        label: keyword,
        kind: getCompletionKind(def),
        detail: def.category ? `[${def.category}]` : undefined,
        documentation: def.description,
        insertTextFormat: InsertTextFormat.Snippet
    };

    // Add snippet with placeholders for arguments
    const effectiveArgCount = Array.isArray(def.argCount) ? def.argCount[0] : (def.argCount ?? 0);
    if (effectiveArgCount > 0) {
        const argCount = effectiveArgCount;
        const placeholders = [];
        for (let i = 0; i < argCount; i++) {
            const argName = def.argTypes?.[i] || `arg${i + 1}`;
            placeholders.push(`\${${i + 1}:${argName}}`);
        }
        item.insertText = `${keyword}\n${placeholders.join(' ')}`;
    } else {
        item.insertText = keyword;
    }

    return item;
}

function getCompletionKind(def: CrystalTagDefinition): CompletionItemKind {
    switch (def.type) {
        case 'block-start':
            return CompletionItemKind.Class;
        case 'keyword':
            return CompletionItemKind.Keyword;
        case 'value':
            return CompletionItemKind.Value;
        default:
            return CompletionItemKind.Text;
    }
}

/**
 * Resolve additional details for a completion item.
 *
 * @param item - Completion item to resolve.
 * @returns Completion item with rich documentation.
 */
export function resolveCrystalCompletion(item: CompletionItem): CompletionItem {
    const def = CRYSTAL_TAGS[item.label];
    if (def) {
        item.documentation = {
            kind: 'markdown',
            value: formatDocumentation(item.label, def)
        };
    }
    return item;
}

function formatDocumentation(keyword: string, def: CrystalTagDefinition): string {
    const lines: string[] = [];

    lines.push(`## ${keyword}`);
    lines.push('');

    if (def.description) {
        lines.push(def.description);
        lines.push('');
    }

    if (def.category) {
        lines.push(`**Category:** ${def.category}`);
    }

    if (def.type) {
        lines.push(`**Type:** ${def.type}`);
    }

    if (def.argCount !== undefined) {
        const count = Array.isArray(def.argCount) ? `${def.argCount[0]}-${def.argCount[1]}` : String(def.argCount);
        lines.push(`**Arguments:** ${count}`);
    }

    if (def.argTypes) {
        lines.push(`**Argument types:** ${def.argTypes.join(', ')}`);
    }

    if (def.defaultValue !== undefined) {
        lines.push(`**Default:** ${def.defaultValue}`);
    }

    if (def.seeAlso && def.seeAlso.length > 0) {
        lines.push(`**See also:** ${def.seeAlso.join(', ')}`);
    }

    return lines.join('\n');
}
