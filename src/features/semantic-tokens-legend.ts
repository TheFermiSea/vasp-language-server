export const tokenTypes = [
    'namespace',
    'type',
    'class',
    'enum',
    'interface',
    'struct',
    'typeParameter',
    'parameter',
    'variable',
    'property',
    'enumMember',
    'event',
    'function',
    'method',
    'macro',
    'keyword',
    'modifier',
    'comment',
    'string',
    'number',
    'regexp',
    'operator'
];

export const tokenModifiers = [
    'declaration',
    'definition',
    'readonly',
    'static',
    'deprecated',
    'abstract',
    'async',
    'modification',
    'documentation',
    'defaultLibrary'
];

export const semanticTokensLegend = {
    tokenTypes,
    tokenModifiers
};

export enum TokenType {
    namespace = 0,
    type = 1,
    class = 2,
    enum = 3,
    interface = 4,
    struct = 5,
    typeParameter = 6,
    parameter = 7,
    variable = 8,
    property = 9,
    enumMember = 10,
    event = 11,
    function = 12,
    method = 13,
    macro = 14,
    keyword = 15,
    modifier = 16,
    comment = 17,
    string = 18,
    number = 19,
    regexp = 20,
    operator = 21
}

export enum TokenModifier {
    declaration = 0,
    definition = 1,
    readonly = 2,
    static = 3,
    deprecated = 4,
    abstract = 5,
    async = 6,
    modification = 7,
    documentation = 8,
    defaultLibrary = 9
}
