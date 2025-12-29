import { SemanticTokensBuilder } from 'vscode-languageserver/node';
import { PoscarDocument } from './parsing';
import { tokenTypeIndices, tokenModifierIndices } from '../semantic-tokens-legend';

export function getPoscarSemanticTokens(parsed: PoscarDocument): { data: number[] } {
    const builder = new SemanticTokensBuilder();

    for (const line of parsed.lines) {
        // Iterate through all tokens in the line
        for (const token of line.tokens) {
            // Mapping Logic
            let typeIdx: number | undefined;
            let modifierIdx = 0; // Default: no modifier

            switch (line.type) {
                case 'comment':
                    typeIdx = tokenTypeIndices.get('comment');
                    break;
                case 'scaling':
                case 'lattice':
                case 'numAtoms':
                case 'positions':
                case 'positionsSelDyn':
                case 'velocities':
                case 'lattVelocitiesVels':
                case 'lattVelocitiesLatt':
                case 'lattVelocitiesState':
                    if (token.type === 'number') {
                        typeIdx = tokenTypeIndices.get('number');
                    } else if (token.type === 'comment') {
                        typeIdx = tokenTypeIndices.get('comment');
                    } else if (token.type === 'constant') {
                        // T/F flags
                        typeIdx = tokenTypeIndices.get('keyword');
                    }
                    break;
                case 'speciesNames':
                    if (token.type === 'string') {
                        typeIdx = tokenTypeIndices.get('class'); // Species are like classes
                        modifierIdx = 1 << (tokenModifierIndices.get('declaration') || 0); // Modifiers are bitflags
                    }
                    break;
                case 'selDynamics':
                case 'positionMode':
                case 'velocityMode':
                case 'lattVelocitiesStart':
                    if (token.type === 'constant') {
                        typeIdx = tokenTypeIndices.get('keyword');
                    } else if (token.type === 'comment') {
                        typeIdx = tokenTypeIndices.get('comment');
                    }
                    break;
            }

            // Fallbacks if specific logic didn't catch it
            if (typeIdx === undefined) {
                if (token.type === 'comment') typeIdx = tokenTypeIndices.get('comment');
                if (token.type === 'number') typeIdx = tokenTypeIndices.get('number');
                if (token.type === 'constant') typeIdx = tokenTypeIndices.get('keyword');
            }

            if (typeIdx !== undefined) {
                builder.push(
                    token.range.start.line,
                    token.range.start.character,
                    token.text.length,
                    typeIdx,
                    modifierIdx
                );
            }
        }
    }

    return builder.build();
}
