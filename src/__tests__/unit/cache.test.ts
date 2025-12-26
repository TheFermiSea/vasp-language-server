import { DocumentCache } from '../../core/document-cache';
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
});
