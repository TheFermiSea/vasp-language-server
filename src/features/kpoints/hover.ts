import { Hover, MarkupKind, Position } from 'vscode-languageserver/node';

/**
 * Provide hover text for KPOINTS files based on line position.
 *
 * @param text - Full document text.
 * @param position - Hover position.
 * @returns Hover payload or null when no hover is available.
 */
export function getKpointsHover(text: string, position: Position): Hover | null {
    const lines = text.split(/\r?\n/);
    if (position.line >= lines.length) return null;

    const lineContent = lines[position.line];
    let contents = '';

    if (position.line === 0) {
        contents = '**KPOINTS Title**\n\nComment line describing the grid.';
    } else if (position.line === 1) {
        contents = '**Number of K-Points**\n\n0 = Automatic Generation.\n>0 = Explicit number (rarely used manually).';
    } else if (position.line === 2) {
        // Mode Line
        if (/^[gG]/.test(lineContent)) {
            contents =
                '**Gamma Centered Grid**\n\nMesh includes the Gamma point (0,0,0). Mandatory for Hexagonal lattices.';
        } else if (/^[mM]/.test(lineContent)) {
            contents =
                '**Monkhorst-Pack Grid**\n\nStandard grid shifted off-Gamma to improve convergence for some cubic systems.';
        } else if (/^[aA]/.test(lineContent)) {
            contents = '**Automatic Mode**\n\nFully automatic generation based on length parameter (R_k).';
        } else if (/^[lL]/.test(lineContent)) {
            contents = '**Line Mode**\n\nGenerates k-points along specific high-symmetry paths for Band Structures.';
        } else {
            contents = '**Generation Mode**\n\nSelects how the grid is built (Gamma, Monkhorst-Pack, Automatic).';
        }
    } else if (position.line >= 3) {
        contents = '**Grid Definition**\n\nDefines the mesh density (Nx Ny Nz) or path coordinates.';
    }

    if (!contents) return null;

    return {
        contents: {
            kind: MarkupKind.Markdown,
            value: contents
        }
    };
}
