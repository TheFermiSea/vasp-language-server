import { TextDocument } from 'vscode-languageserver-textdocument';
import { getIncarSemanticTokens } from '../../features/incar/semantic-tokens';
import { TokenType, TokenModifier } from '../../features/semantic-tokens-legend';

describe('INCAR Semantic Tokens', () => {
    it('should tokenize tags, values, and comments', () => {
        const content = `
ENCUT = 520
# comment
ISMEAR = -5
LREAL = .TRUE.
ALGO = Fast
`;
        const document = TextDocument.create('file:///INCAR', 'vasp', 1, content);
        const tokens = getIncarSemanticTokens(document);

        // We expect a flat array of integers (deltas)
        const data = tokens.data;
        expect(data.length).toBeGreaterThan(0);

        // Decode a few tokens manually to verify logic
        // 1. ENCUT (Line 1, Char 0)
        // DeltaLine: 1, DeltaStart: 0, Len: 5, Type: property (9), Mod: declaration (1)
        expect(data[0]).toBe(1); // Delta Line (from line 0)
        expect(data[1]).toBe(0); // Delta Start
        expect(data[2]).toBe(5); // Length
        expect(data[3]).toBe(TokenType.property);
        expect(data[4]).toBe(1 << TokenModifier.declaration);

        // 2. 520 (Line 1, Char 8) -> DeltaLine: 0, DeltaStart: 8
        expect(data[5]).toBe(0);
        expect(data[6]).toBe(8);
        expect(data[7]).toBe(3); // "520" length
        expect(data[8]).toBe(TokenType.number);
        expect(data[9]).toBe(0);

        // 3. # comment (Line 2)
        // DeltaLine: 1, Type: comment
        // ... find the comment token in the stream
    });

    it('should handle boolean keywords', () => {
        const content = 'LREAL = .TRUE.';
        const document = TextDocument.create('file:///INCAR', 'vasp', 1, content);
        const tokens = getIncarSemanticTokens(document);
        const data = tokens.data;

        // LREAL
        expect(data[3]).toBe(TokenType.property);

        // .TRUE.
        expect(data[8]).toBe(TokenType.keyword);
    });
});
