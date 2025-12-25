import { spawn, ChildProcessWithoutNullStreams } from 'child_process';
import * as path from 'path';

describe('POTCAR Integration', () => {
    let server: ChildProcessWithoutNullStreams;
    let buffer = '';
    const fixtureDir = path.join(__dirname, '../../../test/fixtures/potcar_check');

    beforeEach(() => {
        const serverPath = path.join(__dirname, '../../../out/server.js');
        server = spawn('node', [serverPath, '--stdio']);
        buffer = '';

        server.stdout.on('data', (data) => {
            buffer += data.toString();
        });
    });

    afterEach(async () => {
        if (server && !server.killed) {
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
                    } catch (e) {
                        // ignore parse errors
                    }
                    offset = headerIndex + headerLen + length;
                    if (offset >= buffer.length) break;
                }
            }, 50);

            setTimeout(() => {
                clearInterval(checkInterval);
                reject(new Error('Timeout waiting for message'));
            }, 5000);
        });
    }

    test('Detects mismatch in real file access', async () => {
        // Initialize
        send({
            jsonrpc: "2.0", id: 1, method: "initialize",
            params: { processId: process.pid, rootUri: null, capabilities: {} }
        });
        await waitForMessage((m) => m.id === 1);

        // Initialized
        send({ jsonrpc: "2.0", method: "initialized", params: {} });

        const potcarPath = path.join(fixtureDir, 'POTCAR');

        // Open POTCAR
        send({
            jsonrpc: "2.0",
            method: "textDocument/didOpen",
            params: {
                textDocument: {
                    uri: `file://${potcarPath}`,
                    languageId: "vasp",
                    version: 1,
                    text: "   PAW_PBE O 08Apr2002\n      VRHFIN = O: s2p4\n   End of Potential\n   PAW_PBE Fe 06Sep2000\n      VRHFIN = Fe: d7s1\n"
                }
            }
        });

        // Expect Diagnostics
        const diagMsg = await waitForMessage((m) => m.method === 'textDocument/publishDiagnostics');

        expect(diagMsg.params.uri).toContain("POTCAR");
        const messages = diagMsg.params.diagnostics.map((d: any) => d.message).join(' ');
        // We expect "Mismatch: POSCAR expects 'Fe', but found 'O'"
        expect(messages).toContain("expects 'Fe'");
    });
});
