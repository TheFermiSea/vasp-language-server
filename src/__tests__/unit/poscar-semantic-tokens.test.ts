import { TextDocument } from 'vscode-languageserver-textdocument';
import { parsePoscar } from '../../features/poscar/parsing';
import { getPoscarSemanticTokens } from '../../features/poscar/semantic-tokens';

describe('POSCAR Semantic Tokens', () => {
    test('tokenizes multiple POSCAR sections', () => {
        const content = `Comment
1.0
1 0 0
0 1 0
0 0 1
Fe O
1 1
Selective dynamics
Direct
0 0 0 T F T
Cartesian
0 0 0`;
        const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
        const parsed = parsePoscar(doc);
        const tokens = getPoscarSemanticTokens(parsed);
        expect(tokens.data.length).toBeGreaterThan(0);
    });
});
