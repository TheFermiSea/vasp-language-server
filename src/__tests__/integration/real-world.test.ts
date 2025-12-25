import { spawn, ChildProcessWithoutNullStreams } from 'child_process';
import * as path from 'path';
import * as fs from 'fs';

describe('Real World VASP Files', () => {
    let server: ChildProcessWithoutNullStreams;
    let buffer = '';

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
            await new Promise(r => setTimeout(r, 100));
        }
    });

    function send(msg: any) {
        const str = JSON.stringify(msg);
        const packet = `Content-Length: ${Buffer.byteLength(str, 'utf-8')}\r\n\r\n${str}`;
        server.stdin.write(packet);
    }

    async function waitForMessage(predicate: (msg: any) => boolean): Promise<any> {
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
                        const json = JSON.parse(jsonStr);
                        if (predicate(json)) {
                            clearInterval(checkInterval);
                            resolve(json);
                            return;
                        }
                    } catch (e) { }
                    offset = headerIndex + headerLen + length;
                    if (offset >= buffer.length) break;
                }
            }, 50);

            setTimeout(() => {
                clearInterval(checkInterval);
                reject(new Error('Timeout waiting for message'));
            }, 3000);
        });
    }

    async function validateFile(filePath: string): Promise<any[]> {
        const content = fs.readFileSync(filePath, 'utf-8');
        const uri = `file://${filePath}`;

        // Initialize
        send({ jsonrpc: "2.0", id: 1, method: "initialize", params: { processId: process.pid, rootUri: null, capabilities: {} } });
        await waitForMessage((m) => m.id === 1);
        send({ jsonrpc: "2.0", method: "initialized", params: {} });

        // Open
        send({
            jsonrpc: "2.0",
            method: "textDocument/didOpen",
            params: {
                textDocument: { uri, languageId: "vasp", version: 1, text: content }
            }
        });

        // Wait for Diagnostics
        const msg = await waitForMessage((m) => m.method === 'textDocument/publishDiagnostics' && m.params.uri === uri);
        return msg.params.diagnostics;
    }

    test('Valid: Standard Relaxation (INCAR_RELAX)', async () => {
        const diagnostics = await validateFile(path.join(__dirname, '../../../test/fixtures/real_world/INCAR_RELAX'));
        expect(diagnostics).toHaveLength(0);
    });

    test('Valid: Static Calculation (INCAR_STATIC)', async () => {
        const diagnostics = await validateFile(path.join(__dirname, '../../../test/fixtures/real_world/INCAR_STATIC'));
        expect(diagnostics).toHaveLength(0);
    });

    test('Malformed: Broken Relaxation (INCAR_RELAX_BROKEN)', async () => {
        const diagnostics = await validateFile(path.join(__dirname, '../../../test/fixtures/malformed/INCAR_RELAX_BROKEN'));
        expect(diagnostics.length).toBeGreaterThan(0);

        const messages = diagnostics.map((d: any) => d.message).join(' ');
        expect(messages).toContain("expects a single value"); // ENCUT = 400 eV (parsed as two tokens)
        expect(messages).toContain("Expected integer"); // ISMEAR = Gaussian
        expect(messages).toContain("Expected integer"); // NCORE = Four
    });
});
