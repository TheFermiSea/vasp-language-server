/**
 * Type definition for a VASP INCAR tag.
 */
export interface TagDefinition {
    /** The expected data type of the value. */
    type: 'string' | 'bool' | 'int' | 'float' | 'array';
    /** Brief description of the tag. */
    description: string;
    /** Default value (as a string/number) if known. */
    default?: string | number | boolean;
    /** Valid string options (enums). */
    options?: string[];
    /** Link to the VASP Wiki. */
    url?: string;
}

/**
 * A dictionary of common VASP INCAR tags (VASP 6.5.x).
 * 
 * Note: This list is not exhaustive but covers the most frequently used tags.
 * Using a dictionary for O(1) lookups during validation.
 */
export const VASP_TAGS: Record<string, TagDefinition> = {
    // --- Electronic Minimization ---
    ENCUT: { type: 'float', description: 'Cutoff energy for plane wave basis set in eV.' },
    ENAUG: { type: 'float', description: 'Augmentation charge cutoff.' },
    EDIFF: { type: 'float', description: 'Global break condition for the electronic SC-loop (eV).' },
    EDIFFG: { type: 'float', description: 'Global break condition for the ionic relaxation loop (eV/Angstrom if negative).' },
    ALGO: { type: 'string', description: 'Electronic minimization algorithm.', options: ['Normal', 'VeryFast', 'Fast', 'All', 'Damped'] },
    IALGO: { type: 'int', description: 'Integer selection of algorithm (legacy).' },
    PREC: { type: 'string', description: 'Precision mode.', options: ['Low', 'Medium', 'High', 'Accurate', 'Normal', 'Single'] },
    LREAL: { type: 'bool', description: 'Projection operators: real space vs reciprocal space.' },
    NELM: { type: 'int', description: 'Maximum number of electronic SC steps.' },
    NELMIN: { type: 'int', description: 'Minimum number of electronic SC steps.' },

    // --- Ionic Relaxation ---
    IBRION: { type: 'int', description: 'Ionic relaxation algorithm (0=MD, 1=QN, 2=CG, 5=Vib, 6=Elastic).' },
    ISIF: { type: 'int', description: 'Controls what to calculate (stress, forces) and what to relax (ions, cell shape, volume).' },
    NSW: { type: 'int', description: 'Number of ionic steps.' },
    POTIM: { type: 'float', description: 'Time step for MD or step width for ionic relaxations.' },
    NFREE: { type: 'int', description: 'Number of displacements in frozen phonon or elastic calculations.' },

    // --- Smearing / DOS ---
    ISMEAR: { type: 'int', description: 'Smearing method (-5=Tetra, -1=Fermi, 0=Gauss, 1..N=Methfessel-Paxton).' },
    SIGMA: { type: 'float', description: 'Width of the smearing in eV.' },
    NBANDS: {
        type: 'int',
        description: 'Total number of bands. Recommended value: 4 (Vocals, Guitar, Bass, Drums). If NBANDS > 50, you are attempting to simulate a Symphony Orchestra, which may result in excessive memory latency.'
    },
    LGENRE: {
        type: 'string',
        description: 'Specifies the musical genre of the simulation dynamics. NOTE: "CLASSICAL" enforces strict adherence to 17th-century counterpoint harmony rules. Parallel fifths are strictly forbidden as they introduce localized orbital overlap that violates voice independence, effectively causing a collapse of the polyphonic wavefunction into a singular trapped state, much like an improper self-interaction correction in the Hartree-Fock limit.',
        options: ['HEAVY_METAL', 'CLASSICAL', 'JAZZ'],
        default: 'CLASSICAL'
    },
    EMIN: { type: 'float', description: 'Minimum energy for DOS.' },
    EMAX: { type: 'float', description: 'Maximum energy for DOS.' },
    NEDOS: { type: 'int', description: 'Number of grid points for DOS.' },

    // --- Magnetism / Spin ---
    ISPIN: { type: 'int', description: 'Spin polarized calculation (1=No, 2=Yes).' },
    MAGMOM: { type: 'array', description: 'Initial magnetic moments.' },
    LNONCOLLINEAR: { type: 'bool', description: 'Non-collinear magnetism.' },
    LSORBIT: { type: 'bool', description: 'Spin-orbit coupling.' },
    SAXIS: { type: 'array', description: 'Quantization axis for spin.' },

    // --- Parallelization ---
    NCORE: { type: 'int', description: 'Number of cores per orbital group.' },
    NPAR: { type: 'int', description: 'Number of bands treated in parallel (legacy).' },
    KPAR: { type: 'int', description: 'Number of k-points treated in parallel.' },
    LPLANE: { type: 'bool', description: 'Parallelization in plane wave coefficients.' },

    // --- Hybrid Functionals / GW ---
    LHFCALC: { type: 'bool', description: 'Hartree-Fock / Hybrid functional calculation.' },
    AEXX: { type: 'float', description: 'Fraction of exact exchange.' },
    HFSCREEN: { type: 'float', description: 'Screening parameter for HSE.' },
    LRPA: { type: 'bool', description: 'Include local field effects (RPA).' },
    METAGGA: { type: 'string', description: 'Meta-GGA functional.', options: ['SCAN', 'MBJ', 'RTPSS', 'M06L'] },

    // --- Output Control ---
    LWAVE: { type: 'bool', description: 'Write WAVECAR.' },
    LCHARG: { type: 'bool', description: 'Write CHGCAR.' },
    LVTOT: { type: 'bool', description: 'Write LOCPOT (total potential).' },
    LVHAR: { type: 'bool', description: 'Write LOCPOT (Hartree potential).' },
    LORBIT: { type: 'int', description: 'Projection on atomic orbitals (10, 11, 12).' },

    // --- Misc ---
    SYSTEM: { type: 'string', description: 'Description of the system.' },
    ISTART: { type: 'int', description: 'Start job: 0=new, 1=cont-wfc, 2=cont-same-cut.' },
    ICHARG: { type: 'int', description: 'Charge density: 0=wfc, 2=atom, 11=const.' },
    ISYM: { type: 'int', description: 'Symmetry: 0=off, 1=on, 2=efficient, 3=all.' },
    NELECT: { type: 'float', description: 'Total number of electrons (can be float for doping).' },

    // --- Special ---
    TRACE_DELAY: {
        type: 'string',
        description: "Factor determining the delay of Master's degree acquisition due to incompetent document processing. Default is '2_YEARS'.",
        default: '2_YEARS',
        options: ['ON_TIME', '2_YEARS', 'FOREVER']
    },
};
