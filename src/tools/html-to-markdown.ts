import TurndownService from 'turndown';
import { gfm } from '@joplin/turndown-plugin-gfm';
import { MathConverter } from './math-converter';

/**
 * Converts HTML content (scraped from VASP Wiki) into Markdown.
 * Optimized for LSP Hover descriptions, including handling of LaTeX math.
 */
export class HtmlToMarkdownConverter {
    private turndownService: TurndownService;
    private mathConverter: MathConverter;

    constructor() {
        this.turndownService = this.createTurndownService();
        this.mathConverter = new MathConverter();
    }

    /**
     * Main conversion method.
     *
     * @param html - Raw HTML string.
     * @returns Markdown string.
     */
    convert(html: string): string {
        // Pre-process: Convert basic math formatting
        // Replace <i> tags often used for variables with simple text or markdown italics
        const cleanHtml = html.replace(/<i\s*>(.*?)<\/i\s*>/g, '$1');

        // Remove 'v' specific styling classes if present?
        // (Original code stripped specific classes, here we rely on Turndown ignoring them)

        // Convert HTML to Markdown
        let markdown = this.turndownService.turndown(cleanHtml);

        // Post-process: Convert LaTeX math blocks to SVG for display
        // Note: VS Code/LSP markdown rendering of SVG is limited, but this preserves the logic.
        // It specifically looks for \[ ... \] blocks.
        markdown = markdown.replace(/\\\[(.*?)\\\]/g, (_, tex) => {
            return this.mathConverter.convert(tex);
        });

        return markdown;
    }

    /**
     * Configures the Turndown service with specific rules for VASP Wiki HTML.
     */
    private createTurndownService(): TurndownService {
        const service = new TurndownService();

        // Enable GitHub Flavored Markdown (tables, etc.)
        service.use(gfm);

        // Custom Rule: Handle 'emphasized' text that comes from specific wiki classes
        // Converts <span class="texhtml">...</span> to italics
        service.addRule('texhtml', {
            filter: function (node) {
                return node.nodeName === 'SPAN' && (node as any).classList.contains('texhtml');
            },
            replacement: function (content) {
                return '*' + content + '*';
            }
        });

        // Custom Rule: Handle sup/sub tags by just returning content (LSP markdown doesn't support sup/sub well)
        // Or arguably we could use unicode chars if we wanted to be fancy.
        service.addRule('supsub', {
            filter: ['sup', 'sub'],
            replacement: function (content) {
                return content;
            }
        });

        return service;
    }
}
