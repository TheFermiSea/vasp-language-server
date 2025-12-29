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

/**
 * Caches parsed results (ASTs) for open documents.
 */
export class DocumentCache {
    private cache = new Map<string, { version: number; structure: VaspStructure }>();
    private order: string[] = [];
    private readonly MAX_SIZE = 50;

    public get(document: TextDocument): VaspStructure | undefined {
        const cached = this.cache.get(document.uri);
        if (cached && cached.version === document.version) {
            return cached.structure;
        }
        return undefined;
    }

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
        if (this.order.length > this.MAX_SIZE) {
            const oldest = this.order.shift();
            if (oldest) {
                this.cache.delete(oldest);
            }
        }
    }

    public delete(uri: string): void {
        this.cache.delete(uri);
        const index = this.order.indexOf(uri);
        if (index > -1) {
            this.order.splice(index, 1);
        }
    }
}
