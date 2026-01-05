import { TextDocument } from 'vscode-languageserver-textdocument';
import { Diagnostic, DiagnosticSeverity, Range } from 'vscode-languageserver-types';

export interface PotcarElement {
    symbol: string;
    line: number;
    description: string;
}

export interface ParsedPotcar {
    elements: PotcarElement[];
    diagnostics: Diagnostic[];
}

/**
 * Helper to create a diagnostic with consistent formatting
 */
function createDiagnostic(
    line: number,
    startChar: number,
    endChar: number,
    message: string,
    severity: DiagnosticSeverity
): Diagnostic {
    return {
        severity,
        range: Range.create(line, startChar, line, endChar),
        message
    };
}

/**
 * Known element suffixes in VASP POTCARs that should be stripped to get the base element.
 * Examples: Fe_pv -> Fe, Ti_sv_GW -> Ti, Ge_d -> Ge
 * Common suffixes: _sv (s+p valence), _pv (p valence), _s (soft), _h (hard),
 *                  _d (d states), _GW (GW optimized), _AE (all-electron), _2/_3 (lanthanides)
 */
const ELEMENT_SUFFIXES = ['_sv_GW', '_pv_GW', '_d_GW', '_GW', '_sv', '_pv', '_s', '_h', '_d', '_AE', '_3', '_2'];

/**
 * Extracts the base element symbol from a POTCAR element name by stripping known suffixes.
 * Examples:
 *   "Fe_pv" -> "Fe"
 *   "Ti_sv_GW" -> "Ti"
 *   "O" -> "O"
 *   "Ge_d" -> "Ge"
 *
 * @param elementName - Raw element name from POTCAR metadata.
 * @returns Base element symbol with known suffixes removed.
 */
export function extractBaseElement(elementName: string): string {
    let base = elementName;
    for (const suffix of ELEMENT_SUFFIXES) {
        if (base.endsWith(suffix)) {
            base = base.slice(0, -suffix.length);
            break; // Only strip one suffix (they're ordered from longest to shortest)
        }
    }
    return base;
}

/**
 * Parses a POTCAR file to extract the list of elements.
 * Supports multiple VASP POTCAR format generations:
 *
 * Modern formats (PAW, VASP 5.x+):
 *   - PAW_PBE Fe 06Sep2000
 *   - PAW_LDA Ti_sv 26Sep2005
 *   - PAW_GGA O_s 07Sep2000
 *   - PAW_PBE Ti_sv_GW 05Dec2013
 *   - VRHFIN = Fe: d7s1
 *
 * Legacy formats:
 *   - PAW Fe 06Sep2000 (plain PAW without functional, old LDA)
 *   - US Ti (ultrasoft pseudopotential, minimal header)
 *   - US Fe 08Apr2002 (ultrasoft with date)
 *   - TITEL = PAW_PBE Fe 06Sep2000
 *   - TITEL = US Ti
 *
 * @param document - LSP text document for a POTCAR file.
 * @returns Parsed POTCAR elements and diagnostics.
 */
export function parsePotcar(document: TextDocument): ParsedPotcar {
    const text = document.getText();
    const lines = text.split(/\r?\n/);
    const elements: PotcarElement[] = [];
    const diagnostics: Diagnostic[] = [];

    // Check for empty file
    if (lines.length === 0 || (lines.length === 1 && lines[0].trim() === '')) {
        diagnostics.push(
            createDiagnostic(
                0,
                0,
                0,
                'POTCAR file is empty. A valid POTCAR must contain pseudopotential data with VRHFIN or PAW/US headers identifying element(s).',
                DiagnosticSeverity.Error
            )
        );
        return { elements, diagnostics };
    }

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i];

        // Match VRHFIN = Element : configuration
        // This is the most robust way to find the element identity in modern POTCARs
        // Works for all PAW potentials (PBE, LDA, GGA, GW variants)
        const vrhfinMatch = line.match(/^\s*VRHFIN\s*=\s*([A-Za-z]+)\s*:/);
        if (vrhfinMatch) {
            elements.push({
                symbol: vrhfinMatch[1].trim(),
                line: i,
                description: line.trim()
            });
            continue;
        }
    }

    // If VRHFIN strategy yields nothing (e.g. very old POTCARs or incomplete files), try TITEL/header parsing
    if (elements.length === 0) {
        for (let i = 0; i < lines.length; i++) {
            const line = lines[i];

            // Try multiple formats in order of specificity

            // Format 1: TITEL = PAW_xxx Element Date or TITEL = US Element
            // Example: "   TITEL  = PAW_PBE Fe 06Sep2000" or "   TITEL = US Ti"
            const titelTagMatch = line.match(
                /^\s*TITEL\s*=\s*(?:PAW(?:_(?:PBE|LDA|GGA))?|US)\s+([A-Za-z][A-Za-z0-9_]*)/i
            );
            if (titelTagMatch) {
                const fullName = titelTagMatch[1];
                const baseElement = extractBaseElement(fullName);
                elements.push({
                    symbol: baseElement,
                    line: i,
                    description: line.trim()
                });
                continue;
            }

            // Format 2: Modern PAW with functional: PAW_PBE/PAW_LDA/PAW_GGA Element [Date]
            // Example: "   PAW_PBE Fe 06Sep2000" or "   PAW_PBE Fe_pv 07Sep2000"
            const pawFunctionalMatch = line.match(/^\s*PAW_(?:PBE|LDA|GGA)\s+([A-Za-z][A-Za-z0-9_]*)(?:\s+|$)/);
            if (pawFunctionalMatch) {
                const fullName = pawFunctionalMatch[1];
                const baseElement = extractBaseElement(fullName);
                elements.push({
                    symbol: baseElement,
                    line: i,
                    description: line.trim()
                });
                continue;
            }

            // Format 3: Plain PAW (old LDA without explicit functional): PAW Element [Date]
            // Example: "   PAW Ti_sv 26Sep2005"
            // Must not match "PAW_PBE" etc., so require space after PAW
            const pawPlainMatch = line.match(/^\s*PAW\s+([A-Za-z][A-Za-z0-9_]*)(?:\s+|$)/);
            if (pawPlainMatch && !line.match(/^\s*PAW_/)) {
                const fullName = pawPlainMatch[1];
                const baseElement = extractBaseElement(fullName);
                elements.push({
                    symbol: baseElement,
                    line: i,
                    description: line.trim()
                });
                continue;
            }

            // Format 4: Ultrasoft pseudopotentials (legacy, VASP ~2002): US Element [Date]
            // Example: "   US Ti" or "   US Fe 08Apr2002"
            // These are the oldest format with minimal header info
            const usMatch = line.match(/^\s*US\s+([A-Za-z][A-Za-z0-9_]*)(?:\s+|$)/);
            if (usMatch) {
                const fullName = usMatch[1];
                const baseElement = extractBaseElement(fullName);
                elements.push({
                    symbol: baseElement,
                    line: i,
                    description: line.trim()
                });
                continue;
            }
        }
    }

    // Validation: No elements found
    if (elements.length === 0) {
        diagnostics.push(
            createDiagnostic(
                0,
                0,
                lines[0].length,
                `No valid element definitions found in POTCAR. Expected to find either:\n  - VRHFIN = Element : configuration (modern PAW format)\n  - PAW_PBE/PAW_LDA/PAW_GGA Element Date (PAW header)\n  - US Element (ultrasoft pseudopotential header)\nVerify that this is a valid VASP POTCAR file.`,
                DiagnosticSeverity.Error
            )
        );
    }

    // Validation: Check for common valid element symbols
    const validElements = [
        'H',
        'He',
        'Li',
        'Be',
        'B',
        'C',
        'N',
        'O',
        'F',
        'Ne',
        'Na',
        'Mg',
        'Al',
        'Si',
        'P',
        'S',
        'Cl',
        'Ar',
        'K',
        'Ca',
        'Sc',
        'Ti',
        'V',
        'Cr',
        'Mn',
        'Fe',
        'Co',
        'Ni',
        'Cu',
        'Zn',
        'Ga',
        'Ge',
        'As',
        'Se',
        'Br',
        'Kr',
        'Rb',
        'Sr',
        'Y',
        'Zr',
        'Nb',
        'Mo',
        'Tc',
        'Ru',
        'Rh',
        'Pd',
        'Ag',
        'Cd',
        'In',
        'Sn',
        'Sb',
        'Te',
        'I',
        'Xe',
        'Cs',
        'Ba',
        'La',
        'Ce',
        'Pr',
        'Nd',
        'Pm',
        'Sm',
        'Eu',
        'Gd',
        'Tb',
        'Dy',
        'Ho',
        'Er',
        'Tm',
        'Yb',
        'Lu',
        'Hf',
        'Ta',
        'W',
        'Re',
        'Os',
        'Ir',
        'Pt',
        'Au',
        'Hg',
        'Tl',
        'Pb',
        'Bi',
        'Po',
        'At',
        'Rn',
        'Fr',
        'Ra',
        'Ac',
        'Th',
        'Pa',
        'U',
        'Np',
        'Pu',
        'Am',
        'Cm',
        'Bk',
        'Cf',
        'Es',
        'Fm',
        'Md',
        'No',
        'Lr',
        'Rf',
        'Db',
        'Sg',
        'Bh',
        'Hs',
        'Mt',
        'Ds',
        'Rg',
        'Cn',
        'Nh',
        'Fl',
        'Mc',
        'Lv',
        'Ts',
        'Og'
    ];

    for (const element of elements) {
        if (!validElements.includes(element.symbol)) {
            diagnostics.push(
                createDiagnostic(
                    element.line,
                    0,
                    lines[element.line].length,
                    `Line ${element.line + 1}: Unrecognized element symbol '${element.symbol}'. This may be a parsing error or an unusual element name. Standard element symbols are 1-2 letters (e.g., 'Fe', 'O', 'Ti').`,
                    DiagnosticSeverity.Warning
                )
            );
        }
    }

    // Validation: Check for duplicate elements (warning, as this might be intentional)
    const elementCounts = new Map<string, number>();
    for (const element of elements) {
        elementCounts.set(element.symbol, (elementCounts.get(element.symbol) || 0) + 1);
    }
    for (const [symbol, count] of elementCounts) {
        if (count > 1) {
            diagnostics.push(
                createDiagnostic(
                    0,
                    0,
                    lines[0].length,
                    `Element '${symbol}' appears ${count} times in POTCAR. This is unusual - typically each element appears once. Verify that the POTCAR matches your POSCAR species order.`,
                    DiagnosticSeverity.Warning
                )
            );
        }
    }

    return { elements, diagnostics };
}
