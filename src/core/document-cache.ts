import { TextDocument } from 'vscode-languageserver-textdocument';
import { IncarDocument } from '../features/incar/parsing';
import { PoscarDocument } from '../features/poscar/parsing';
import { KpointsDocument } from '../features/kpoints/parsing';
import { CrystalDocument } from '../features/crystal/parsing';

/**
 * Union of all supported DFT document structures.
 * Supports VASP (INCAR, POSCAR, KPOINTS, POTCAR) and CRYSTAL23 (.d12)
 */
export type DFTStructure =
    | { type: 'incar'; data: IncarDocument }
    | { type: 'poscar'; data: PoscarDocument }
    | { type: 'kpoints'; data: KpointsDocument }
    | { type: 'potcar'; data: null }
    | { type: 'crystal'; data: CrystalDocument };

// Backwards compatibility alias
export type VaspStructure = DFTStructure;

/** Default maximum number of documents to cache */
export const DEFAULT_CACHE_SIZE = 50;

/**
 * Caches parsed results (ASTs) for open documents.
 */
export class DocumentCache {
    private cache = new Map<string, { version: number; structure: VaspStructure }>();
    private order: string[] = [];
    private readonly maxSize: number;

    /**
     * Creates a new DocumentCache.
     * @param maxSize Maximum number of documents to cache (default: 50)
     */
    constructor(maxSize: number = DEFAULT_CACHE_SIZE) {
        this.maxSize = maxSize;
    }

    /**
     * Retrieve a cached structure for a document if the version matches.
     *
     * @param document - The document to look up in the cache.
     * @returns The cached structure when present and up-to-date, otherwise undefined.
     */
    public get(document: TextDocument): VaspStructure | undefined {
        const cached = this.cache.get(document.uri);
        if (cached && cached.version === document.version) {
            return cached.structure;
        }
        return undefined;
    }

    /**
     * Store a parsed structure in the cache and enforce FIFO eviction.
     *
     * @param document - The document providing the cache key and version.
     * @param structure - Parsed document structure to cache.
     */
    public set(document: TextDocument, structure: VaspStructure): void {
        const uri = document.uri;

        if (!this.cache.has(uri)) {
            this.order.push(uri);
        }

        this.cache.set(uri, {
            version: document.version,
            structure
        });

        // Enforce FIFO limit
        if (this.order.length > this.maxSize) {
            const oldest = this.order.shift();
            if (oldest) {
                this.cache.delete(oldest);
            }
        }
    }

    /**
     * Remove a document from the cache and eviction queue.
     *
     * @param uri - Document URI to remove.
     */
    public delete(uri: string): void {
        this.cache.delete(uri);
        const index = this.order.indexOf(uri);
        if (index > -1) {
            this.order.splice(index, 1);
        }
    }
}
