import { TextDocument } from 'vscode-languageserver-textdocument';

export interface PotcarElement {
    symbol: string;
    line: number;
    description: string;
}

export interface ParsedPotcar {
    elements: PotcarElement[];
}

/**
 * Parses a POTCAR file to extract the list of elements.
 * Looks for the "VRHFIN" tag or the first line of each potential block.
 * VASP POTCAR headers often look like:
 *    PAW_PBE Fe 06Sep2000
 *    VRHFIN = Fe: d7s1
 */
export function parsePotcar(document: TextDocument): ParsedPotcar {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
    const elements: PotcarElement[] = [];

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i];

        // Match VRHFIN = Element : configuration
        // This is a robust way to find the element identity in modern POTCARs
        const vrhfinMatch = line.match(/^\s*VRHFIN\s*=\s*([A-Za-z]+)\s*:/);
        if (vrhfinMatch) {
            elements.push({
                symbol: vrhfinMatch[1].trim(),
                line: i,
                description: line.trim()
            });
            continue;
        }

        // Also check for TITEL fields (PAW_PBE Element Date)
        // This usually appears at the very top of a block
        // Regex: space PAW_... Element Date
        const titelMatch = line.match(/^\s*PAW_(?:PBE|LDA|GGA)\s+([A-Za-z0-9_]+)\s+/);
        if (titelMatch) {
            // We prioritize VRHFIN if we found one recently? 
            // Actually, a POTCAR block has both. We should avoid duplicates.
            // But VRHFIN is inside the block. TITEL is at start. 
            // Let's stick to VRHFIN as it's definitive for the element properties.
            // OR, if VASP 5.x POTKUL, TITEL is good.
        }
    }

    // If VRHFIN strategy yields nothing (e.g. very old POTCARs), try TITEL
    // If VRHFIN strategy yields nothing (e.g. very old POTCARs), try TITEL
    if (elements.length === 0) {
        for (let i = 0; i < lines.length; i++) {
            const line = lines[i];
            // Example: "   PAW_PBE Fe 06Sep2000"
            // Careful not to match random text. 
            // TITEL in POTCAR usually starts with spaces.
            const titleMatch = line.match(/^\s*PAW_(?:PBE|LDA|GGA)\s+([A-Za-z0-9]+)(?:\s+|$)/);
            if (titleMatch) {
                // If the element name has suffix like Fe_pv, extract just Fe or keep Fe_pv?
                // Usually POSCAR just says "Fe".
                // Let's keep strict match for now, user can see mismatch if it occurs.
                elements.push({
                    symbol: titleMatch[1],
                    line: i,
                    description: lines[i].trim()
                });
            }
        }
    }

    return { elements };
}
