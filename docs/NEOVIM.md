# Neovim Setup Guide
>
> **Recommended**: Use [LazyVim](https://lazyvim.org) or `lazy.nvim`.

## ⚡ Accelerated Setup

### 1. Install the Server

```bash
npm install -g vasp-language-server
```

### 2. Install the Neovim Plugin (Recommended)

Using `lazy.nvim`:

```lua
return {
  {
    "TheFermiSea/vasp-language-server",
    ft = "vasp",
    config = function()
      require("vasp-lsp").setup()
    end,
  },
}
```

### 2. Configure LazyVim

Add this to `lua/plugins/vasp.lua`:

```lua
return {
  {
    "neovim/nvim-lspconfig",
    opts = {
      servers = {
        vasp_ls = {
          -- Defaults to looking for 'vasp-lsp' in PATH
          cmd = { "vasp-lsp", "--stdio" }, 
          filetypes = { "vasp" },
        },
      },
      setup = {
        vasp_ls = function()
          -- Add filetype detection for VASP files
          vim.filetype.add({
            pattern = {
              ['.*INCAR.*'] = 'vasp',
              ['.*POSCAR.*'] = 'vasp',
              ['.*CONTCAR.*'] = 'vasp',
              ['.*POTCAR.*'] = 'vasp',
              ['.*KPOINTS.*'] = 'vasp',
            },
          })
        end,
      },
    },
  },
}
```

---

## 🔌 Recommended Plugins

Maximize your productivity with these plugins:

1. **[nvim-cmp](https://github.com/hrsh7th/nvim-cmp)**: For intelligent autocompletion of VASP tags.
2. **[telescope.nvim](https://github.com/nvim-telescope/telescope.nvim)**: Use `:Telescope lsp_document_symbols` to navigate large INCAR/POSCAR files.
3. **[catppuccin](https://github.com/catppuccin/nvim)**: A theme that supports our Semantic Tokens (colors booleans, numbers, and tags distinctly).

---

## 📦 Releasing the Neovim Plugin

See [NEOVIM_RELEASE.md](NEOVIM_RELEASE.md) for the release checklist.

---

## 🛠 Manual Setup (Vanilla `init.lua`)

If you don't use a plugin manager:

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

-- 2. Setup LSP
local lspconfig = require('lspconfig')
local configs = require('lspconfig.configs')

if not configs.vasp_ls then
  configs.vasp_ls = {
    default_config = {
      cmd = { 'vasp-lsp', '--stdio' },
      filetypes = { 'vasp' },
      root_dir = lspconfig.util.root_pattern('.git', 'INCAR', 'POSCAR'),
      settings = {},
    },
  }
end

lspconfig.vasp_ls.setup{}
```
