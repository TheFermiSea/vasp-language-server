import { Connection } from 'vscode-languageserver';

export class Logger {
    private connection: Connection | null = null;

    initialize(connection: Connection): void {
        this.connection = connection;
    }

    info(message: string): void {
        this.log('INFO', message);
    }

    warn(message: string): void {
        this.log('WARN', message);
    }

    error(message: string, error?: unknown): void {
        const errorMessage = error instanceof Error ? `: ${error.message}` : '';
        this.log('ERROR', `${message}${errorMessage}`);
        if (error instanceof Error && error.stack) {
            this.log('ERROR', error.stack);
        }
    }

    private log(level: string, message: string): void {
        if (!this.connection) {
            // Fallback to console if not initialized (though these might not show up anywhere useful)
            console.log(`[${level}] ${message}`);
            return;
        }
        const timestamp = new Date().toISOString();
        // LSP protocol uses connection.console.log/.error/.warn
        if (level === 'ERROR') {
            this.connection.console.error(`[${timestamp}] ${message}`);
        } else if (level === 'WARN') {
            this.connection.console.warn(`[${timestamp}] ${message}`);
        } else {
            this.connection.console.log(`[${timestamp}] ${message}`);
        }
    }
}

export const logger = new Logger();
