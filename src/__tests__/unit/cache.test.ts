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

    test('should enforce FIFO limit of 50', () => {
        const docs: TextDocument[] = [];
        for (let i = 0; i < 60; i++) {
            docs.push(TextDocument.create(`file:///test${i}.incar`, 'vasp', 1, ''));
            cache.set(docs[i], { type: 'incar' as const, data: { statements: [], allTokens: [] } });
        }

        // First 10 should be gone (60 - 50 = 10 removed)
        for (let i = 0; i < 10; i++) {
            expect(cache.get(docs[i])).toBeUndefined();
        }

        // Last 50 should remain
        for (let i = 10; i < 60; i++) {
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

        // Add 7 documents to a cache with size 5
        for (let i = 0; i < 7; i++) {
            docs.push(TextDocument.create(`file:///custom${i}.incar`, 'vasp', 1, ''));
            customCache.set(docs[i], { type: 'incar' as const, data: { statements: [], allTokens: [] } });
        }

        // First 2 should be evicted (7 - 5 = 2 removed)
        expect(customCache.get(docs[0])).toBeUndefined();
        expect(customCache.get(docs[1])).toBeUndefined();

        // Last 5 should remain
        for (let i = 2; i < 7; i++) {
            expect(customCache.get(docs[i])).toBeDefined();
        }
    });

    test('should use default size when no parameter provided', () => {
        const defaultCache = new DocumentCache();
        const docs: TextDocument[] = [];

        // Add 55 documents - should evict first 5
        for (let i = 0; i < 55; i++) {
            docs.push(TextDocument.create(`file:///default${i}.incar`, 'vasp', 1, ''));
            defaultCache.set(docs[i], { type: 'incar' as const, data: { statements: [], allTokens: [] } });
        }

        // First 5 should be evicted (55 - 50 = 5 removed)
        for (let i = 0; i < 5; i++) {
            expect(defaultCache.get(docs[i])).toBeUndefined();
        }

        // Last 50 should remain
        for (let i = 5; i < 55; i++) {
            expect(defaultCache.get(docs[i])).toBeDefined();
        }
    });
});
