import { Mwn } from 'mwn';
import {
    extractCategoryPageIds,
    fetchIncarTags,
    parseIncarTagsFromPages,
    PagesResponse,
    CategoryMembersResponse
} from '../../tools/vasp-wiki';

jest.mock('mwn', () => ({
    Mwn: jest.fn()
}));

describe('VASP Wiki Scraper', () => {
    const mwnMock = Mwn as unknown as jest.Mock;
    const converter = { convert: (html: string) => html };

    afterEach(() => {
        jest.resetAllMocks();
    });

    test('extractCategoryPageIds returns page ids', () => {
        const response: CategoryMembersResponse = {
            query: { categorymembers: [{ pageid: 1 }, { pageid: 2 }] }
        };
        expect(extractCategoryPageIds(response)).toEqual([1, 2]);
    });

    test('parseIncarTagsFromPages skips missing fields', () => {
        const pages = {
            a: { title: 'ENCUT' },
            b: { text: { '*': '<p>Missing title</p>' } }
        };
        const tags = parseIncarTagsFromPages(pages, converter);
        expect(tags).toHaveLength(0);
    });

    test('parseIncarTagsFromPages builds tags from HTML', () => {
        const pages = {
            a: { title: 'ENCUT', text: { '*': '<p>Energy cutoff</p>' } }
        };
        const tags = parseIncarTagsFromPages(pages, converter);
        expect(tags).toHaveLength(1);
        expect(tags[0].name).toBe('ENCUT');
        expect(tags[0].description).toContain('Energy cutoff');
    });

    test('fetchIncarTags returns empty when category members missing', async () => {
        const requestMock = jest.fn().mockResolvedValue({ query: {} });
        mwnMock.mockImplementation(() => ({ request: requestMock }));

        const tags = await fetchIncarTags('https://example.com');
        expect(tags).toEqual([]);
        expect(requestMock).toHaveBeenCalledTimes(1);
    });

    test('fetchIncarTags returns empty when pages missing', async () => {
        const requestMock = jest
            .fn()
            .mockResolvedValueOnce({ query: { categorymembers: [{ pageid: 1 }] } })
            .mockResolvedValueOnce({ query: {} });
        mwnMock.mockImplementation(() => ({ request: requestMock }));

        const tags = await fetchIncarTags('https://example.com');
        expect(tags).toEqual([]);
        expect(requestMock).toHaveBeenCalledTimes(2);
    });

    test('fetchIncarTags returns tags on success', async () => {
        const pagesResponse: PagesResponse = {
            query: { pages: { a: { title: 'ENCUT', text: { '*': '<p>Energy cutoff</p>' } } } }
        };
        const requestMock = jest
            .fn()
            .mockResolvedValueOnce({ query: { categorymembers: [{ pageid: 1 }] } })
            .mockResolvedValueOnce(pagesResponse);
        mwnMock.mockImplementation(() => ({ request: requestMock }));

        const tags = await fetchIncarTags('https://example.com');
        expect(tags).toHaveLength(1);
        expect(tags[0].name).toBe('ENCUT');
    });

    test('fetchIncarTags handles request errors', async () => {
        const requestMock = jest.fn().mockRejectedValue(new Error('Timeout'));
        mwnMock.mockImplementation(() => ({ request: requestMock }));

        const tags = await fetchIncarTags('https://example.com');
        expect(tags).toEqual([]);
        expect(requestMock).toHaveBeenCalledTimes(1);
    });
});
