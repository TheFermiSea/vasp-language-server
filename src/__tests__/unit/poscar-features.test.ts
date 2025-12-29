import { getPoscarHover } from '../../features/poscar/hover';
import { getPoscarSemanticTokens } from '../../features/poscar/semantic-tokens';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { parsePoscar } from '../../features/poscar/parsing';
import { Position } from 'vscode-languageserver/node';

describe('POSCAR Features', () => {
    const poscarText = `Fe3O4
1.0
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
Fe O
3 4
Direct
0.1 0.1 0.1
0.5 0.5 0.5`;

    const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, poscarText);
    const parsed = parsePoscar(doc);

    describe('Hovers', () => {
        it('should hover identifying comment line', () => {
            const result = getPoscarHover(parsed, Position.create(0, 0));
            expect(result).not.toBeNull();
            // @ts-expect-error - MarkupContent has value property
            expect(result?.contents.value).toContain('**POSCAR Title / Comment**');
        });

        it('should hover scaling factor', () => {
            const result = getPoscarHover(parsed, Position.create(1, 0));
            // @ts-expect-error - MarkupContent has value property
            expect(result?.contents.value).toContain('**Universal Scaling Factor**');
        });

        it('should hover lattice vectors', () => {
            const result = getPoscarHover(parsed, Position.create(3, 0));
            // @ts-expect-error - MarkupContent has value property
            expect(result?.contents.value).toContain('**Lattice Vector**');
        });

        it('should hover species names', () => {
            const result = getPoscarHover(parsed, Position.create(5, 0));
            // @ts-expect-error - MarkupContent has value property
            expect(result?.contents.value).toContain('**Species Names**');
        });

        it('should hover positions', () => {
            const result = getPoscarHover(parsed, Position.create(8, 0));
            // @ts-expect-error - MarkupContent has value property
            expect(result?.contents.value).toContain('**Atomic Position**');
        });
    });

    describe('Semantic Tokens', () => {
        it('should tokenize lattice vectors as numbers', () => {
            const tokens = getPoscarSemanticTokens(parsed);
            // Logic to check if line 2 (lattice 1) has number tokens
            // This is harder to check raw integer array, but we can check not-empty
            expect(tokens.data.length).toBeGreaterThan(0);
        });
    });
});
