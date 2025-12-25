import { CodeAction, CodeActionKind, Diagnostic, TextEdit } from "vscode-languageserver/node";
import { VASP_TAGS } from "./data/vasp-tags";
import { levenshteinDistance } from "./util";

/**
 * Generates Code Actions (Quick Fixes) for INCAR diagnostics.
 * 
 * @param uri - The document URI.
 * @param diagnostics - The diagnostics in the current context.
 * @returns An array of CodeAction items.
 */
export function getIncarCodeActions(uri: string, diagnostics: Diagnostic[]): CodeAction[] {
    const actions: CodeAction[] = [];

    for (const diagnostic of diagnostics) {
        // Look for "Unknown tag" messages
        if (diagnostic.message.startsWith("Unknown tag")) {
            // Extract the unknown tag name from the message
            // Pattern: "Unknown tag 'TAGNAME'. Check spelling..."
            const match = diagnostic.message.match(/Unknown tag '(.+)'/);
            if (match && match[1]) {
                const unknownTag = match[1];

                // Find closest match in VASP_TAGS
                let closestTag = '';
                let minDist = Infinity;

                for (const validTag of Object.keys(VASP_TAGS)) {
                    const dist = levenshteinDistance(unknownTag, validTag);
                    if (dist < minDist) {
                        minDist = dist;
                        closestTag = validTag;
                    }
                }

                // Threshold for suggestion (e.g. distance <= 3)
                // Also ensures we don't suggest if the distance is too large relative to length
                if (minDist <= 3 && closestTag) {
                    const action: CodeAction = {
                        title: `Change to '${closestTag}'`,
                        kind: CodeActionKind.QuickFix,
                        diagnostics: [diagnostic],
                        edit: {
                            changes: {
                                [uri]: [
                                    TextEdit.replace(diagnostic.range, closestTag)
                                ]
                            }
                        }
                    };
                    actions.push(action);
                }
            }
        }
    }

    return actions;
}
