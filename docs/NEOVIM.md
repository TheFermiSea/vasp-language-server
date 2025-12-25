
# Neovim Setup Guide for VASP Language Server

The VASP Language Server is built on the standard Language Server Protocol (LSP), making it fully compatible with Neovim's built-in LSP client (requires Neovim 0.5+, recommended 0.9+ for Semantic Tokens).

## 1. Prerequisites

- **Node.js**: The server runs on Node.js.
- **npm**: To install the server.

## 2. Installation

While we work on adding this to `mason.nvim`, you can install it manually:

```bash
# From the source directory
npm install
npm run compile
# The executable is now at ./out/server.js
```

## 3. Configuration (lua)

Add the following to your `init.lua` or `plugins/lsp.lua`. This configuration handles:

1. **Filetype Detection**: Automatically detecting `INCAR`, `POSCAR`, `POTCAR`, `KPOINTS`.
2. **Server Setup**: Launching the server via `stdio`.

```lua
-- 1. Register Filetypes
vim.filetype.add({
  pattern = {
    ['.*INCAR.*'] = 'vasp',
    ['.*POSCAR.*'] = 'vasp',
    ['.*CONTCAR.*'] = 'vasp',
    ['.*POTCAR.*'] = 'vasp',
    ['.*KPOINTS.*'] = 'vasp',
  },
})

-- 2. Configure LSP
local lspconfig = require('lspconfig')
local configs = require('lspconfig.configs')

if not configs.vasp_ls then
  configs.vasp_ls = {
    default_config = {
      cmd = { 'node', '/path/to/vasp-language-server/out/server.js', '--stdio' },
      filetypes = { 'vasp' },
      root_dir = lspconfig.util.root_pattern('.git', 'INCAR', 'POSCAR'),
      settings = {},
    },
  }
end

lspconfig.vasp_ls.setup({
  on_attach = function(client, bufnr)
    -- Enable Semantic Tokens (Highlighter)
    client.server_capabilities.semanticTokensProvider = {
      full = true,
      legend = {
        tokenTypes = { "property", "keyword", "number", "comment", "string" },
        tokenModifiers = { "declaration" }
      },
      range = false,
    }
  end,
})
```

## 4. Feature Checklist

| Feature | Neovim Support | Notes |
| :--- | :--- | :--- |
| **Diagnostics** | ✅ Built-in | Works out of the box with `vim.diagnostic`. |
| **Auto-Completion** | ✅ Built-in | Works with `nvim-cmp` or `omnifunc`. |
| **Hover Docs** | ✅ Built-in | Use `K` or `vim.lsp.buf.hover()`. |
| **Code Actions** | ✅ Built-in | Use `vim.lsp.buf.code_action()` (fixes typos). |
| **Outline** | ✅ Built-in | Use `Telescope lsp_document_symbols` or `Aerial`. |
| **Semantic Tokens** | ✅ Built-in (0.9+) | Requires a theme supporting semantic highlights (e.g., Catppuccin, Tokyo Night). |
| **Folding** | ✅ Built-in | Use `set foldmethod=expr foldexpr=nvim_treesitter#foldexpr()` or LSP folding capabilities. |
