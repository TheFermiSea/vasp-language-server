import { Hover, MarkupKind, Position } from 'vscode-languageserver/node';
import { PoscarDocument, PoscarBlockType } from './parsing';

export function getPoscarHover(parsed: PoscarDocument, position: Position): Hover | null {
    // Find the line corresponding to the hover position
    const line = parsed.lines.find((l) => l.line.lineNumber === position.line);
    if (!line) return null;

    let contents = '';

    switch (line.type) {
        case 'comment':
            contents = '**POSCAR Title / Comment**\n\nIdentifying string for the structure.';
            break;
        case 'scaling':
            contents = '**Universal Scaling Factor**\n\nMultiplies all lattice vectors. Usually 1.0.';
            break;
        case 'lattice':
            contents = '**Lattice Vector**\n\nBasis vector definition in Ångstroms.';
            break;
        case 'speciesNames':
            contents = '**Species Names**\n\nOrdered list of atomic species (e.g., Fe O). Must match POTCAR.';
            break;
        case 'numAtoms':
            contents = '**Ion Counts**\n\nNumber of atoms for each species defined above.';
            break;
        case 'selDynamics':
            contents = '**Selective Dynamics Mode**\n\nEnables T/F flags for constraining atom movement.';
            break;
        case 'positionMode':
            contents = '**Coordinate System**\n\nDefines if positions are Direct (fractional) or Cartesian (Å).';
            break;
        case 'positions':
        case 'positionsSelDyn':
            contents = '**Atomic Position**\n\nCoordinates of an ion.';
            break;
        case 'velocities':
            contents = '**Initial Velocity**\n\nInitial velocity for molecular dynamics.';
            break;
    }

    if (!contents) return null;

    return {
        contents: {
            kind: MarkupKind.Markdown,
            value: contents
        }
    };
}
