import { LspServer } from './core/lsp-server';

/**
 * Entry point for the VASP Language Server.
 */
const server = new LspServer();
server.start();
