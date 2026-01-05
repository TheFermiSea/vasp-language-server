import * as vscode from 'vscode';
import { LanguageClient, LanguageClientOptions, ServerOptions } from 'vscode-languageclient/node';

let client: LanguageClient | undefined;

export function activate(context: vscode.ExtensionContext): void {
    const config = vscode.workspace.getConfiguration('vaspLsp');
    const serverPath = config.get<string>('serverPath', 'vasp-lsp');

    const serverOptions: ServerOptions = {
        command: serverPath,
        args: ['--stdio']
    };

    const clientOptions: LanguageClientOptions = {
        documentSelector: [
            { scheme: 'file', language: 'vasp' },
            { scheme: 'untitled', language: 'vasp' }
        ],
        outputChannelName: 'VASP Language Server'
    };

    client = new LanguageClient('vaspLanguageServer', 'VASP Language Server', serverOptions, clientOptions);
    context.subscriptions.push(client.start());
}

export async function deactivate(): Promise<void> {
    if (client) {
        await client.stop();
        client = undefined;
    }
}
