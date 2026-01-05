import { TextDocument } from 'vscode-languageserver-textdocument';
import { Position } from 'vscode-languageserver/node';
import { parsePoscar } from '../../features/poscar/parsing';
import { getPoscarHover } from '../../features/poscar/hover';

describe('POSCAR Hover', () => {
    const content = `Title
1.0
1 0 0
0 1 0
0 0 1
Fe
1
Selective dynamics
Direct
0 0 0 T T F
Cartesian
0 0 0`;

    const doc = TextDocument.create('file:///POSCAR', 'vasp', 1, content);
    const parsed = parsePoscar(doc);

    test('hovers selective dynamics line', () => {
        const hover = getPoscarHover(parsed, Position.create(7, 0));
        expect(hover?.contents).toBeDefined();
        const value = (hover?.contents as { value: string }).value;
        expect(value).toContain('Selective Dynamics');
    });

    test('hovers position mode line', () => {
        const hover = getPoscarHover(parsed, Position.create(8, 0));
        expect(hover?.contents).toBeDefined();
        const value = (hover?.contents as { value: string }).value;
        expect(value).toContain('Coordinate System');
    });

    test('hovers positions with selective dynamics flags', () => {
        const hover = getPoscarHover(parsed, Position.create(9, 0));
        expect(hover?.contents).toBeDefined();
        const value = (hover?.contents as { value: string }).value;
        expect(value).toContain('Atomic Position');
    });

    test('hovers velocities line', () => {
        const hover = getPoscarHover(parsed, Position.create(11, 0));
        expect(hover?.contents).toBeDefined();
        const value = (hover?.contents as { value: string }).value;
        expect(value).toContain('Initial Velocity');
    });
});
