/** @type {import('ts-jest').JestConfigWithTsJest} */
module.exports = {
    preset: 'ts-jest',
    testEnvironment: 'node',
    roots: ['<rootDir>/src'],
    testMatch: [
        "**/__tests__/**/*.+(ts|tsx|js)",
        "**/?(*.)+(spec|test).+(ts|tsx|js)"
    ],
    transform: {
        "^.+\\.(ts|tsx)$": "ts-jest"
    },
    moduleFileExtensions: ['ts', 'tsx', 'js', 'jsx', 'json', 'node'],
    verbose: true,
    // Helper to allow absolute imports if needed, though relative is fine for now
    moduleDirectories: ['node_modules', 'src'],
    // Coverage thresholds to maintain code quality
    // These are based on tested source files only (not entry points/tools)
    coverageThreshold: {
        global: {
            statements: 80,
            branches: 65,
            functions: 80,
            lines: 80
        }
    },
    // Exclude non-essential files from coverage
    coveragePathIgnorePatterns: [
        '/node_modules/',
        '/src/__tests__/',
        '/src/tools/',
        '/src/server.ts'
    ]
};
