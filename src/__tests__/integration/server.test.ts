import { spawn, ChildProcessWithoutNullStreams } from 'child_process';
import * as path from 'path';

describe('LSP Server Integration', () => {
    let server: ChildProcessWithoutNullStreams;
    let buffer = '';

    type JsonRpcMessage = {
        jsonrpc?: string;
        id?: number | string;
        method?: string;
        params?: Record<string, unknown>;
        result?: unknown;
        error?: unknown;
    };

    type PublishDiagnosticsParams = {
        uri: string;
        diagnostics: Array<{ message: string; range: Record<string, unknown> }>;
    };

    beforeEach(() => {
        const serverPath = path.join(__dirname, '../../../out/server.js');
        server = spawn('node', [serverPath, '--stdio']);
        buffer = '';

        server.stdout.on('data', (data) => {
            buffer += data.toString();
        });
    });

    afterEach(async () => {
        if (!server.killed) {
            server.stdin.end();
            server.kill();
            // Give it a moment to die to avoid race conditions with next test
            await new Promise((r) => setTimeout(r, 100));
        }
    });

    function send(msg: JsonRpcMessage) {
        const str = JSON.stringify(msg);
        const packet = `Content-Length: ${Buffer.byteLength(str, 'utf-8')}\r\n\r\n${str}`;
        server.stdin.write(packet);
    }

    async function waitForMessage(predicate: (msg: JsonRpcMessage) => boolean): Promise<JsonRpcMessage> {
        return new Promise((resolve, reject) => {
            const checkInterval = setInterval(() => {
                let offset = 0;
                while (true) {
                    const match = buffer.slice(offset).match(/Content-Length: (\d+)\r\n\r\n/);
                    if (!match) break;

                    const length = parseInt(match[1]);
                    const headerLen = match[0].length;
                    const headerIndex = buffer.indexOf(match[0], offset);

                    if (buffer.length < headerIndex + headerLen + length) break;

                    const jsonStr = buffer.slice(headerIndex + headerLen, headerIndex + headerLen + length);
                    try {
                        const json = JSON.parse(jsonStr) as JsonRpcMessage;
                        if (predicate(json)) {
                            clearInterval(checkInterval);
                            clearTimeout(timeoutId);
                            resolve(json);
                            return;
                        }
                    } catch {
                        // ignore parse errors
                    }
                    offset = headerIndex + headerLen + length;
                    if (offset >= buffer.length) break;
                }
            }, 50);

            const timeoutId = setTimeout(() => {
                clearInterval(checkInterval);
                reject(new Error('Timeout waiting for message'));
            }, 3000);
        });
    }

    test('initializes and returns capabilities', async () => {
        send({
            jsonrpc: '2.0',
            id: 1,
            method: 'initialize',
            params: { processId: process.pid, rootUri: null, capabilities: {} }
        });

        const msg = await waitForMessage((m) => m.id === 1);
        const result = msg.result as { capabilities?: { textDocumentSync?: unknown; hoverProvider?: boolean } };
        expect(result.capabilities?.textDocumentSync).toBeDefined();
        expect(result.capabilities?.hoverProvider).toBe(true);
    });

    test('validates INCAR file on open', async () => {
        // Initialize first
        send({
            jsonrpc: '2.0',
            id: 1,
            method: 'initialize',
            params: { processId: process.pid, rootUri: null, capabilities: {} }
        });
        await waitForMessage((m) => m.id === 1);

        // Send Initialized
        send({ jsonrpc: '2.0', method: 'initialized', params: {} });

        // Open Document with Error
        send({
            jsonrpc: '2.0',
            method: 'textDocument/didOpen',
            params: {
                textDocument: {
                    uri: 'file:///test/INCAR',
                    languageId: 'vasp',
                    version: 1,
                    text: 'ENCUT = High'
                }
            }
        });

        // Expect Diagnostics
        const diagMsg = await waitForMessage((m) => m.method === 'textDocument/publishDiagnostics');
        const diagParams = diagMsg.params as PublishDiagnosticsParams;
        expect(diagParams.uri).toBe('file:///test/INCAR');
        expect(diagParams.diagnostics).toHaveLength(1);
        expect(diagParams.diagnostics[0].message).toContain('Expected number');
    });

    test('provides Document Symbols (Outline) for INCAR', async () => {
        // Init
        send({
            jsonrpc: '2.0',
            id: 1,
            method: 'initialize',
            params: { processId: process.pid, rootUri: null, capabilities: {} }
        });
        await waitForMessage((m) => m.id === 1);
        send({ jsonrpc: '2.0', method: 'initialized', params: {} });

        // Open
        send({
            jsonrpc: '2.0',
            method: 'textDocument/didOpen',
            params: {
                textDocument: {
                    uri: 'file:///test/INCAR_SYMBOLS',
                    languageId: 'vasp',
                    version: 1,
                    text: 'ENCUT = 500\nISMEAR = 0'
                }
            }
        });

        // Wait for document to be processed (diagnostics published)
        await waitForMessage(
            (m) =>
                m.method === 'textDocument/publishDiagnostics' &&
                (m.params as PublishDiagnosticsParams | undefined)?.uri === 'file:///test/INCAR_SYMBOLS'
        );

        // Request Symbols
        send({
            jsonrpc: '2.0',
            id: 2,
            method: 'textDocument/documentSymbol',
            params: { textDocument: { uri: 'file:///test/INCAR_SYMBOLS' } }
        });

        const msg = await waitForMessage((m) => m.id === 2);
        if (msg.error) {
            console.error('LSP Error:', JSON.stringify(msg.error, null, 2));
        }
        const symbols = msg.result as Array<{ name: string }>;
        expect(symbols).toBeDefined();
        expect(symbols).toHaveLength(2);
        expect(symbols[0].name).toBe('ENCUT');
        expect(symbols[1].name).toBe('ISMEAR');
    });

    test('provides Code Actions for typos', async () => {
        // Init
        send({
            jsonrpc: '2.0',
            id: 1,
            method: 'initialize',
            params: { processId: process.pid, rootUri: null, capabilities: {} }
        });
        await waitForMessage((m) => m.id === 1);
        send({ jsonrpc: '2.0', method: 'initialized', params: {} });

        const docUri = 'file:///test/INCAR_TYPO';
        // Open
        send({
            jsonrpc: '2.0',
            method: 'textDocument/didOpen',
            params: {
                textDocument: { uri: docUri, languageId: 'vasp', version: 1, text: 'EMCUT = 500' }
            }
        });

        // Wait for Diagnostic
        const diagMsg = await waitForMessage(
            (m) =>
                m.method === 'textDocument/publishDiagnostics' &&
                (m.params as PublishDiagnosticsParams | undefined)?.uri === docUri
        );
        const diagParams = diagMsg.params as PublishDiagnosticsParams;
        const diagnostic = diagParams.diagnostics[0];
        expect(diagnostic.message).toContain("Unknown tag 'EMCUT'");

        // Request CodeAction
        send({
            jsonrpc: '2.0',
            id: 2,
            method: 'textDocument/codeAction',
            params: {
                textDocument: { uri: docUri },
                range: diagnostic.range,
                context: { diagnostics: [diagnostic] }
            }
        });

        const actionMsg = await waitForMessage((m) => m.id === 2);
        const actions = actionMsg.result as Array<{ title: string; kind?: string }>;
        expect(actions.length).toBeGreaterThan(0);
        expect(actions[0].title).toBe("Change to 'ENCUT'");
        expect(actions[0].kind).toBe('quickfix');
    });
});
