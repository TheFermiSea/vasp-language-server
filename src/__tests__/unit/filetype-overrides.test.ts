import { LspServer } from '../../core/lsp-server';

function createFakeConnection() {
    return {
        console: {
            log: jest.fn(),
            warn: jest.fn(),
            error: jest.fn()
        },
        onInitialize: jest.fn(),
        onInitialized: jest.fn(),
        onCompletion: jest.fn(),
        onCompletionResolve: jest.fn(),
        onHover: jest.fn(),
        onDocumentFormatting: jest.fn(),
        onCodeAction: jest.fn(),
        onDocumentSymbol: jest.fn(),
        onFoldingRanges: jest.fn(),
        onExecuteCommand: jest.fn(),
        onDidChangeConfiguration: jest.fn(),
        onDidOpenTextDocument: jest.fn(),
        onDidChangeTextDocument: jest.fn(),
        onDidCloseTextDocument: jest.fn(),
        onDidSaveTextDocument: jest.fn(),
        onWillSaveTextDocument: jest.fn(),
        onWillSaveTextDocumentWaitUntil: jest.fn(),
        sendDiagnostics: jest.fn(),
        workspace: {
            getConfiguration: jest.fn().mockResolvedValue({})
        },
        languages: {
            semanticTokens: {
                on: jest.fn()
            }
        },
        listen: jest.fn()
    };
}

jest.mock('vscode-languageserver/node', () => {
    const actual = jest.requireActual('vscode-languageserver/node');
    return {
        ...actual,
        createConnection: jest.fn(() => createFakeConnection())
    };
});

describe('File type overrides', () => {
    test('uses filename and extension overrides', () => {
        const server = new LspServer({
            fileTypeOverrides: {
                filenames: {
                    'INCAR.relax': 'incar'
                },
                extensions: {
                    '.vasp': 'poscar',
                    kp: 'kpoints'
                }
            }
        });

        const getFileType = (server as unknown as { getFileType: (uri: string) => string }).getFileType.bind(
            server
        );

        expect(getFileType('file:///tmp/INCAR.relax')).toBe('incar');
        expect(getFileType('file:///tmp/structure.vasp')).toBe('poscar');
        expect(getFileType('file:///tmp/mesh.kp')).toBe('kpoints');
    });
});
