import { CompletionItem, CompletionItemKind } from 'vscode-languageserver/node';
import { VASP_TAGS } from '../../data/vasp-tags';

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

export function resolveIncarCompletion(item: CompletionItem): CompletionItem {
    const tagDef = VASP_TAGS[item.data as string];
    if (tagDef) {
        item.documentation = tagDef.description + (tagDef.default ? `\n\nDefault: ${tagDef.default}` : '');
    }
    return item;
}
