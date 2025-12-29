/**
 * CRYSTAL23 keyword definitions
 *
 * This file contains definitions for CRYSTAL23 input file keywords,
 * including their types, argument counts, and documentation.
 *
 * Categories:
 * - structure: Geometry type and structure definitions
 * - geometry: Geometry optimization keywords
 * - scf: SCF control keywords
 * - dft: DFT functional and exchange-correlation keywords
 * - basis: Basis set related keywords
 * - properties: Post-SCF property calculations
 * - output: Output control keywords
 */

export interface CrystalTagDefinition {
    type: 'keyword' | 'block-start' | 'block-end' | 'value';
    description: string;
    category?: 'structure' | 'geometry' | 'scf' | 'dft' | 'basis' | 'properties' | 'output';
    argCount?: number | [number, number]; // Fixed or [min, max]
    argTypes?: string[]; // Type hints for arguments
    defaultValue?: unknown;
    deprecated?: boolean;
    seeAlso?: string[];
}

export const CRYSTAL_TAGS: Record<string, CrystalTagDefinition> = {
    // ========== Structure Keywords ==========
    CRYSTAL: {
        type: 'keyword',
        category: 'structure',
        description: '3D periodic system with full crystalline symmetry (space groups 1-230)'
    },
    SLAB: {
        type: 'keyword',
        category: 'structure',
        description: '2D periodic system (surface/thin film) with layer groups 1-80'
    },
    POLYMER: {
        type: 'keyword',
        category: 'structure',
        description: '1D periodic system (chain/wire) with rod groups'
    },
    MOLECULE: {
        type: 'keyword',
        category: 'structure',
        description: 'Non-periodic molecular system in Cartesian coordinates'
    },
    HELIX: {
        type: 'keyword',
        category: 'structure',
        description: 'Helical system (nanotubes, DNA, etc.)'
    },
    EXTERNAL: {
        type: 'keyword',
        category: 'structure',
        description: 'Read geometry from external file (fort.34 / .gui format)'
    },
    SUPERCEL: {
        type: 'keyword',
        category: 'structure',
        description: 'Create supercell from primitive cell',
        argCount: [3, 9],
        argTypes: ['int', 'int', 'int']
    },
    NANOTUBE: {
        type: 'keyword',
        category: 'structure',
        description: 'Define nanotube from (n,m) indices',
        argCount: 2,
        argTypes: ['int', 'int']
    },

    // ========== SCF Keywords ==========
    SHRINK: {
        type: 'keyword',
        category: 'scf',
        description: 'K-point mesh: IS (shrinking factor for Pack-Monkhorst net) and IP (Gilat net for Fermi energy)',
        argCount: 2,
        argTypes: ['int', 'int'],
        seeAlso: ['TOLINTEG', 'FMIXING']
    },
    TOLINTEG: {
        type: 'keyword',
        category: 'scf',
        description: 'Integral tolerance thresholds. Five integers controlling Coulomb/exchange overlap screening',
        argCount: 5,
        argTypes: ['int', 'int', 'int', 'int', 'int'],
        defaultValue: [6, 6, 6, 6, 12]
    },
    FMIXING: {
        type: 'keyword',
        category: 'scf',
        description: 'Fock/KS matrix mixing percentage (0-100) for SCF convergence',
        argCount: 1,
        argTypes: ['int'],
        defaultValue: 30
    },
    BROYDEN: {
        type: 'keyword',
        category: 'scf',
        description: 'Broyden mixing for SCF acceleration. Arguments: W0, IMIX, [ISTART]',
        argCount: [2, 3],
        argTypes: ['float', 'int', 'int']
    },
    DIIS: {
        type: 'keyword',
        category: 'scf',
        description: 'Direct Inversion of Iterative Subspace for SCF acceleration',
        argCount: 0
    },
    ANDERSON: {
        type: 'keyword',
        category: 'scf',
        description: 'Anderson mixing for SCF convergence',
        argCount: 0
    },
    MAXCYCLE: {
        type: 'keyword',
        category: 'scf',
        description: 'Maximum number of SCF cycles',
        argCount: 1,
        argTypes: ['int'],
        defaultValue: 50
    },
    LEVSHIFT: {
        type: 'keyword',
        category: 'scf',
        description: 'Level shift to improve SCF convergence. Arguments: SHIFT, ILOCK',
        argCount: 2,
        argTypes: ['float', 'int']
    },
    SPINLOCK: {
        type: 'keyword',
        category: 'scf',
        description: 'Lock spin configuration during SCF',
        argCount: 2,
        argTypes: ['int', 'int']
    },
    TOLDEE: {
        type: 'keyword',
        category: 'scf',
        description: 'Energy convergence threshold (10^-N Hartree)',
        argCount: 1,
        argTypes: ['int'],
        defaultValue: 6
    },
    TOLDEP: {
        type: 'keyword',
        category: 'scf',
        description: 'Density matrix convergence threshold',
        argCount: 1,
        argTypes: ['int']
    },
    TOLDEG: {
        type: 'keyword',
        category: 'scf',
        description: 'Gradient convergence threshold for geometry optimization',
        argCount: 1,
        argTypes: ['float']
    },
    ENDSCF: {
        type: 'block-end',
        category: 'scf',
        description: 'End of SCF block'
    },

    // ========== DFT Keywords ==========
    DFT: {
        type: 'block-start',
        category: 'dft',
        description: 'Begin DFT block. Specify exchange-correlation functional',
        seeAlso: ['B3LYP', 'PBE', 'HSE06']
    },
    ENDDFT: {
        type: 'block-end',
        category: 'dft',
        description: 'End of DFT block'
    },
    B3LYP: {
        type: 'keyword',
        category: 'dft',
        description: 'B3LYP hybrid functional (20% HF exchange)'
    },
    PBE: {
        type: 'keyword',
        category: 'dft',
        description: 'PBE GGA functional (Perdew-Burke-Ernzerhof)'
    },
    PBESOL: {
        type: 'keyword',
        category: 'dft',
        description: 'PBEsol functional (revised PBE for solids)'
    },
    HSE06: {
        type: 'keyword',
        category: 'dft',
        description: 'HSE06 screened hybrid functional. Recommended for band gaps',
        seeAlso: ['HFEXCHANGE']
    },
    M06: {
        type: 'keyword',
        category: 'dft',
        description: 'M06 meta-GGA hybrid functional'
    },
    PW91: {
        type: 'keyword',
        category: 'dft',
        description: 'PW91 GGA functional'
    },
    LDA: {
        type: 'keyword',
        category: 'dft',
        description: 'Local Density Approximation (Dirac-Slater exchange)'
    },
    LSDA: {
        type: 'keyword',
        category: 'dft',
        description: 'Local Spin Density Approximation'
    },
    SPIN: {
        type: 'keyword',
        category: 'dft',
        description: 'Enable spin-polarized calculation',
        argCount: 0
    },
    HFEXCHANGE: {
        type: 'keyword',
        category: 'dft',
        description: 'Fraction of Hartree-Fock exchange in hybrid functional',
        argCount: 1,
        argTypes: ['float']
    },
    XLGRID: {
        type: 'keyword',
        category: 'dft',
        description: 'Extra-large integration grid for DFT'
    },
    XXLGRID: {
        type: 'keyword',
        category: 'dft',
        description: 'Extra-extra-large integration grid for high precision DFT'
    },
    EXCHANGE: {
        type: 'keyword',
        category: 'dft',
        description: 'Specify exchange functional'
    },
    CORRELAT: {
        type: 'keyword',
        category: 'dft',
        description: 'Specify correlation functional'
    },

    // ========== Geometry Optimization ==========
    OPTGEOM: {
        type: 'block-start',
        category: 'geometry',
        description: 'Begin geometry optimization block',
        seeAlso: ['FULLOPTG', 'ATOMONLY', 'CELLONLY']
    },
    ENDOPT: {
        type: 'block-end',
        category: 'geometry',
        description: 'End of geometry optimization block'
    },
    FULLOPTG: {
        type: 'keyword',
        category: 'geometry',
        description: 'Full geometry optimization (cell parameters + atomic positions)',
        seeAlso: ['ATOMONLY', 'CELLONLY']
    },
    ATOMONLY: {
        type: 'keyword',
        category: 'geometry',
        description: 'Optimize atomic positions only, fixed cell',
        seeAlso: ['FULLOPTG', 'CELLONLY']
    },
    CELLONLY: {
        type: 'keyword',
        category: 'geometry',
        description: 'Optimize cell parameters only, fixed atomic fractional coordinates',
        seeAlso: ['FULLOPTG', 'ATOMONLY']
    },
    CVOLOPT: {
        type: 'keyword',
        category: 'geometry',
        description: 'Optimize at constant volume'
    },
    INTRANS: {
        type: 'keyword',
        category: 'geometry',
        description: 'Transform from internal to Cartesian coordinates'
    },
    HESSOPT: {
        type: 'keyword',
        category: 'geometry',
        description: 'Read/write Hessian for geometry optimization restart'
    },
    MAXOPT: {
        type: 'keyword',
        category: 'geometry',
        description: 'Maximum geometry optimization steps',
        argCount: 1,
        argTypes: ['int'],
        defaultValue: 50
    },

    // ========== Frequency & Elastic ==========
    FREQCALC: {
        type: 'block-start',
        category: 'properties',
        description: 'Frequency calculation (phonons/vibrations)'
    },
    ELASTCON: {
        type: 'block-start',
        category: 'properties',
        description: 'Elastic constants calculation'
    },
    EOS: {
        type: 'block-start',
        category: 'properties',
        description: 'Equation of state calculation'
    },

    // ========== Properties ==========
    BAND: {
        type: 'block-start',
        category: 'properties',
        description: 'Band structure calculation along k-path'
    },
    DOSS: {
        type: 'block-start',
        category: 'properties',
        description: 'Density of states calculation'
    },
    ECHG: {
        type: 'keyword',
        category: 'properties',
        description: 'Print electron charge density'
    },
    POTC: {
        type: 'keyword',
        category: 'properties',
        description: 'Electrostatic potential calculation'
    },
    COOP: {
        type: 'keyword',
        category: 'properties',
        description: 'Crystal Orbital Overlap Population analysis'
    },
    BWDF: {
        type: 'keyword',
        category: 'properties',
        description: 'Bloch wave decomposition'
    },
    ANISOTRO: {
        type: 'keyword',
        category: 'properties',
        description: 'Anisotropic temperature factors'
    },

    // ========== Basis Set ==========
    BASISSET: {
        type: 'block-start',
        category: 'basis',
        description: 'Begin basis set definition block'
    },
    GHOSTS: {
        type: 'keyword',
        category: 'basis',
        description: 'Define ghost atoms (basis functions without nuclear charge)',
        argCount: 1,
        argTypes: ['int']
    },
    ECP: {
        type: 'keyword',
        category: 'basis',
        description: 'Effective Core Potential definition'
    },

    // ========== Output Control ==========
    NOSYMADA: {
        type: 'keyword',
        category: 'output',
        description: 'Disable symmetry adaptation of orbitals'
    },
    NOCOMM: {
        type: 'keyword',
        category: 'output',
        description: 'Suppress communication output'
    },
    SETPRINT: {
        type: 'keyword',
        category: 'output',
        description: 'Set print options',
        argCount: 1,
        argTypes: ['int']
    },

    // ========== Termination ==========
    END: {
        type: 'block-end',
        description: 'End of section or entire input'
    },
    ENDG: {
        type: 'block-end',
        category: 'structure',
        description: 'End of geometry section'
    },
    ENDB: {
        type: 'block-end',
        category: 'basis',
        description: 'End of basis set section'
    }
};

/**
 * Get all keywords for a specific category
 */
export function getKeywordsByCategory(category: CrystalTagDefinition['category']): string[] {
    return Object.entries(CRYSTAL_TAGS)
        .filter(([_, def]) => def.category === category)
        .map(([keyword, _]) => keyword);
}

/**
 * Check if a keyword exists (case-insensitive)
 */
export function isValidKeyword(keyword: string): boolean {
    return keyword.toUpperCase() in CRYSTAL_TAGS;
}
