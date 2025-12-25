import { MarkupContent, MarkupKind } from 'vscode-languageserver-types';

/**
 * Represents a single VASP INCAR tag (e.g., "ENCUT", "PREC").
 * Stores its documentation and provides methods to format it for LSP responses.
 */
export class IncarTag {
    /** The tag name (e.g., "ENCUT"). */
    public readonly name: string;

    /** The raw markdown description. */
    public readonly description: string;

    /**
     * @param name - The name of the tag.
     * @param description - The full markdown description.
     */
    constructor(name: string, description: string) {
        this.name = name;
        this.description = description;
    }

    /**
     * Factory method to create an IncarTag from a markdown string.
     * 
     * @param markdown - The markdown content describing the tag.
     * @param tagName - The name of the tag.
     * @returns A new IncarTag instance.
     */
    static fromMarkdown(markdown: string, tagName: string): IncarTag {
        return new IncarTag(tagName, markdown);
    }

    /**
     * Generates a Hover object (MarkupContent) for the LSP.
     * 
     * @param baseUrl - (Optional) Base URL to prepend to relative links in the documentation.
     * @returns A MarkupContent object with kind 'markdown'.
     */
    getHoverText(baseUrl: string): MarkupContent {
        // We could prepend the baseUrl to links here if logic dictates
        // For now, return the stored description
        return {
            kind: MarkupKind.Markdown,
            value: `**${this.name}**\n\n${this.description}`
        };
    }
}
