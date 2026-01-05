import { Range } from 'vscode-languageserver-types';

export interface BaseToken {
    text: string;
    range: Range;
}

export type IncarTokenType = 'tag' | 'equals' | 'value' | 'comment' | 'semicolon' | 'eol' | 'invalid';

export interface IncarToken extends BaseToken {
    type: IncarTokenType;
}

export const poscarTokenTypes = ['comment', 'string', 'number', 'constant', 'invalid'] as const;

export type PoscarTokenType = (typeof poscarTokenTypes)[number];

export interface PoscarToken extends BaseToken {
    type?: PoscarTokenType;
}

export type CrystalTokenType =
    | 'title'
    | 'geometry-type'
    | 'space-group'
    | 'lattice-param'
    | 'atom-count'
    | 'atom-line'
    | 'keyword'
    | 'block-start'
    | 'block-end'
    | 'basis-header'
    | 'basis-data'
    | 'number'
    | 'comment'
    | 'end'
    | 'unknown';

export interface CrystalToken extends BaseToken {
    type: CrystalTokenType;
    lineNumber: number;
}

export type TokenType = IncarTokenType | PoscarTokenType | CrystalTokenType;
