import { CompletionItem, CompletionItemKind } from 'vscode-languageserver/node';
import { VASP_TAGS } from '../../data/vasp-tags';

/**
 * Build completion items for all known INCAR tags.
 *
 * @returns Completion items derived from the VASP tag database.
 */
export function getIncarCompletions(): CompletionItem[] {
    const items: CompletionItem[] = [];
    for (const [tag, def] of Object.entries(VASP_TAGS)) {
        items.push({
            label: tag,
            kind: CompletionItemKind.Property,
            data: tag,
            detail: def.description
        });
    }
    return items;
}

/**
 * Enrich a selected INCAR completion item with documentation details.
 *
 * @param item - Completion item to resolve.
 * @returns Completion item with documentation populated when available.
 */
export function resolveIncarCompletion(item: CompletionItem): CompletionItem {
    const tagDef = VASP_TAGS[item.data as string];
    if (tagDef) {
        item.documentation = tagDef.description + (tagDef.default ? `\n\nDefault: ${tagDef.default}` : '');
    }
    return item;
}
