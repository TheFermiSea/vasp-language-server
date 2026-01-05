import { Position } from 'vscode-languageserver/node';
import { getKpointsHover } from '../../features/kpoints/hover';

describe('KPOINTS Hover', () => {
    test('hovers Gamma mode', () => {
        const content = `Comment\n0\nGamma\n4 4 4\n0 0 0`;
        const hover = getKpointsHover(content, Position.create(2, 0));
        const value = (hover?.contents as { value: string }).value;
        expect(value).toContain('Gamma Centered');
    });

    test('hovers Monkhorst-Pack mode', () => {
        const content = `Comment\n0\nMonkhorst-Pack\n4 4 4\n0 0 0`;
        const hover = getKpointsHover(content, Position.create(2, 0));
        const value = (hover?.contents as { value: string }).value;
        expect(value).toContain('Monkhorst-Pack');
    });

    test('hovers Automatic mode', () => {
        const content = `Comment\n0\nAuto\n20`;
        const hover = getKpointsHover(content, Position.create(2, 0));
        const value = (hover?.contents as { value: string }).value;
        expect(value).toContain('Automatic Mode');
    });

    test('hovers Line mode', () => {
        const content = `Comment\n0\nLine-mode\n0 0 0`;
        const hover = getKpointsHover(content, Position.create(2, 0));
        const value = (hover?.contents as { value: string }).value;
        expect(value).toContain('Line Mode');
    });

    test('hovers default mode hint', () => {
        const content = `Comment\n0\nZ\n4 4 4`;
        const hover = getKpointsHover(content, Position.create(2, 0));
        const value = (hover?.contents as { value: string }).value;
        expect(value).toContain('Generation Mode');
    });

    test('hovers grid definition line', () => {
        const content = `Comment\n0\nGamma\n4 4 4\n0 0 0`;
        const hover = getKpointsHover(content, Position.create(3, 0));
        const value = (hover?.contents as { value: string }).value;
        expect(value).toContain('Grid Definition');
    });
});
