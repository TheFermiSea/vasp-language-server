import { getKpointsCompletions } from '../../features/kpoints/completion';
import { getKpointsHover } from '../../features/kpoints/hover';
import { getKpointsSemanticTokens } from '../../features/kpoints/semantic-tokens';
import { Position } from 'vscode-languageserver/node';

describe('KPOINTS Features', () => {
    const kpointsText = `Comment
0
Gamma
4 4 1
0 0 0`;

    describe('Completion', () => {
        it('should suggestion modes on line 2', () => {
            const items = getKpointsCompletions(kpointsText, Position.create(2, 0));
            expect(items).toHaveLength(6);
            const labels = items.map((i) => i.label);
            expect(labels).toContain('Gamma');
            expect(labels).toContain('Monkhorst-Pack');
            expect(labels).toContain('Snippet: Monkhorst-Pack 4x4x4');
        });

        it('should NOT suggest on other lines', () => {
            const items = getKpointsCompletions(kpointsText, Position.create(0, 0));
            expect(items).toHaveLength(0);
        });
    });

    describe('Hover', () => {
        it('should explain Gamma mode', () => {
            const result = getKpointsHover(kpointsText, Position.create(2, 0));
            // @ts-expect-error - contents.value check
            expect(result?.contents.value).toContain('**Gamma Centered Grid**');
        });

        it('should explain grid definition', () => {
            const result = getKpointsHover(kpointsText, Position.create(3, 0));
            // @ts-expect-error
            expect(result?.contents.value).toContain('**Grid Definition**');
        });
    });

    describe('Semantic Tokens', () => {
        it('should return tokens', () => {
            const tokens = getKpointsSemanticTokens(kpointsText);
            expect(tokens.data.length).toBeGreaterThan(0);
        });
    });
});
