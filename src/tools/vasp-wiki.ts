import { Mwn } from 'mwn';
import { HtmlToMarkdownConverter } from './html-to-markdown';
import { IncarTag } from '../features/incar/tag';
import { logger } from '../utils/logger';

/**
 * Fetches and parses INCAR tag documentation from the VASP Wiki using the MediaWiki API.
 *
 * @param baseUrl - The base URL of the wiki (e.g., https://www.vasp.at/wiki/index.php).
 * @returns A list of parsed IncarTag objects.
 */
export async function fetchIncarTags(baseUrl: string): Promise<IncarTag[]> {
    logger.info('Connecting to VASP Wiki...');

    // Initialize MediaWiki client
    const bot = new Mwn({
        apiUrl: `${baseUrl}/api.php`,
        userAgent: 'VASP-Language-Server/0.1 (https://github.com/yourusername/vasp-lsp)'
    });

    const converter = new HtmlToMarkdownConverter();
    const incarTags: IncarTag[] = [];

    try {
        // Fetch the list of pages in the 'INCAR_tags' category
        const categoryMembers = await bot.request({
            action: 'query',
            list: 'categorymembers',
            cmtitle: 'Category:INCAR_tags',
            cmlimit: 5000
        });

        if (!categoryMembers?.query?.categorymembers) {
            logger.warn('No category members found for INCAR_tags');
            return [];
        }

        // Extract page IDs
        const pageIds = categoryMembers.query.categorymembers.map((member: any) => member.pageid);

        // Fetch the actual content (text) of these pages
        const pagesData = await bot.request({
            action: 'query',
            pageids: pageIds.join('|'),
            prop: 'text'
        });

        if (!pagesData?.query?.pages) {
            logger.warn('No pages data returned');
            return [];
        }

        const pages = pagesData.query.pages;

        // Iterate through each page response
        for (const pageId in pages) {
            const page = pages[pageId];
            const title = page.title;
            const htmlContent = page.text['*'];

            // For now, convert the raw HTML to Markdown.
            const markdown = converter.convert(htmlContent);

            // Create the IncarTag object
            const tag = IncarTag.fromMarkdown(markdown, title);
            incarTags.push(tag);
        }

        logger.info(`Fetched ${incarTags.length} INCAR tags from VASP Wiki.`);
    } catch (error) {
        logger.error('Failed to fetch INCAR tags from wiki', error);
    }

    return incarTags;
}

if (require.main === module) {
    (async () => {
        const tags = await fetchIncarTags('https://www.vasp.at/wiki/index.php');
        if (tags.length > 0) {
            console.log(
                `Fetched ${tags.length} tags. (Saving implementation skipped to avoid overwriting handmade list with raw data for now)`
            );
            // In a real scenario, we would generate the TS file here.
            // For this task, we just wanted to prove the integration exists.
        }
    })();
}
