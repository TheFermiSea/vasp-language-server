import { Connection } from 'vscode-languageserver';
import { logger } from '../../utils/logger';

describe('Logger', () => {
    const resetLogger = () => {
        (logger as unknown as { connection: Connection | null }).connection = null;
    };

    afterEach(() => {
        resetLogger();
        jest.restoreAllMocks();
    });

    test('logs to console when not initialized', () => {
        const spy = jest.spyOn(console, 'log').mockImplementation(() => {});
        logger.info('hello');
        expect(spy).toHaveBeenCalled();
        expect(spy.mock.calls[0][0]).toContain('[INFO] hello');
    });

    test('logs info/warn/error to connection console when initialized', () => {
        const connection = {
            console: {
                log: jest.fn(),
                warn: jest.fn(),
                error: jest.fn()
            }
        } as unknown as Connection;

        logger.initialize(connection);
        logger.info('info');
        logger.warn('warn');
        logger.error('error');

        expect(connection.console.log).toHaveBeenCalled();
        expect(connection.console.warn).toHaveBeenCalled();
        expect(connection.console.error).toHaveBeenCalled();
    });

    test('logs error stack when provided', () => {
        const errorMock = jest.fn();
        const connection = {
            console: {
                log: jest.fn(),
                warn: jest.fn(),
                error: errorMock
            }
        } as unknown as Connection;

        logger.initialize(connection);
        const err = new Error('boom');
        err.stack = 'STACK';
        logger.error('failure', err);

        expect(errorMock).toHaveBeenCalledTimes(2);
        const secondCall = errorMock.mock.calls[1][0] as string;
        expect(secondCall).toContain('STACK');
    });
});
