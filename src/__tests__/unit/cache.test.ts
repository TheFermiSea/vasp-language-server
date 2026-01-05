import { DocumentCache, DEFAULT_CACHE_SIZE } from '../../core/document-cache';
import { TextDocument } from 'vscode-languageserver-textdocument';

describe('DocumentCache', () => {
    let cache: DocumentCache;

    beforeEach(() => {
        cache = new DocumentCache();
    });

    test('should store and retrieve a document', () => {
        const doc = TextDocument.create('file:///test.incar', 'vasp', 1, 'ENCUT = 500');
        const structure = { type: 'incar' as const, data: { statements: [], allTokens: [], errors: [] } };

        cache.set(doc, structure);
        expect(cache.get(doc)).toBe(structure);
    });

    test('should return undefined for version mismatch', () => {
        const doc1 = TextDocument.create('file:///test.incar', 'vasp', 1, 'ENCUT = 500');
        const doc2 = TextDocument.create('file:///test.incar', 'vasp', 2, 'ENCUT = 500');
        const structure = { type: 'incar' as const, data: { statements: [], allTokens: [], errors: [] } };

        cache.set(doc1, structure);
        expect(cache.get(doc2)).toBeUndefined();
    });

    test('should enforce LRU limit of 50', () => {
        const docs: TextDocument[] = [];
        for (let i = 0; i < 50; i++) {
            docs.push(TextDocument.create(`file:///test${i}.incar`, 'vasp', 1, ''));
            cache.set(docs[i], { type: 'incar' as const, data: { statements: [], allTokens: [] } });
        }

        // Touch the first document so it becomes most recently used.
        expect(cache.get(docs[0])).toBeDefined();

        for (let i = 50; i < 60; i++) {
            docs.push(TextDocument.create(`file:///test${i}.incar`, 'vasp', 1, ''));
            cache.set(docs[i], { type: 'incar' as const, data: { statements: [], allTokens: [] } });
        }

        // Documents 1-10 should be gone (least recently used)
        for (let i = 1; i <= 10; i++) {
            expect(cache.get(docs[i])).toBeUndefined();
        }

        // The touched doc0 should remain, along with the most recent additions.
        expect(cache.get(docs[0])).toBeDefined();
        for (let i = 11; i < 60; i++) {
            expect(cache.get(docs[i])).toBeDefined();
        }
    });

    test('should delete entry correctly', () => {
        const doc = TextDocument.create('file:///test.incar', 'vasp', 1, '');
        cache.set(doc, { type: 'incar' as const, data: { statements: [], allTokens: [] } });

        cache.delete(doc.uri);
        expect(cache.get(doc)).toBeUndefined();
    });

    test('should export default cache size constant', () => {
        expect(DEFAULT_CACHE_SIZE).toBe(50);
    });

    test('should support custom cache size', () => {
        const customCache = new DocumentCache(5);
        const docs: TextDocument[] = [];

        // Add 5 documents to a cache with size 5
        for (let i = 0; i < 5; i++) {
            docs.push(TextDocument.create(`file:///custom${i}.incar`, 'vasp', 1, ''));
            customCache.set(docs[i], { type: 'incar' as const, data: { statements: [], allTokens: [] } });
        }

        // Touch first entry to mark it as most recently used.
        expect(customCache.get(docs[0])).toBeDefined();

        // Add 2 more documents to trigger eviction
        for (let i = 5; i < 7; i++) {
            docs.push(TextDocument.create(`file:///custom${i}.incar`, 'vasp', 1, ''));
            customCache.set(docs[i], { type: 'incar' as const, data: { statements: [], allTokens: [] } });
        }

        // Oldest untouched entries should be evicted
        expect(customCache.get(docs[1])).toBeUndefined();
        expect(customCache.get(docs[2])).toBeUndefined();

        // Most recent entries should remain
        expect(customCache.get(docs[0])).toBeDefined();
        for (let i = 3; i < 7; i++) {
            expect(customCache.get(docs[i])).toBeDefined();
        }
    });

    test('should use default size when no parameter provided', () => {
        const defaultCache = new DocumentCache();
        const docs: TextDocument[] = [];

        // Add 50 documents (fill the cache)
        for (let i = 0; i < 50; i++) {
            docs.push(TextDocument.create(`file:///default${i}.incar`, 'vasp', 1, ''));
            defaultCache.set(docs[i], { type: 'incar' as const, data: { statements: [], allTokens: [] } });
        }

        // Touch an early entry so it survives eviction.
        expect(defaultCache.get(docs[0])).toBeDefined();

        // Add 5 more documents - should evict least recently used
        for (let i = 50; i < 55; i++) {
            docs.push(TextDocument.create(`file:///default${i}.incar`, 'vasp', 1, ''));
            defaultCache.set(docs[i], { type: 'incar' as const, data: { statements: [], allTokens: [] } });
        }

        // Entries 1-5 should be evicted (least recently used)
        for (let i = 1; i <= 5; i++) {
            expect(defaultCache.get(docs[i])).toBeUndefined();
        }

        // Touched entry and new ones should remain
        expect(defaultCache.get(docs[0])).toBeDefined();
        for (let i = 6; i < 55; i++) {
            expect(defaultCache.get(docs[i])).toBeDefined();
        }
    });
});
