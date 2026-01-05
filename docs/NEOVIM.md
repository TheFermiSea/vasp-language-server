# Neovim Setup Guide
>
> **Recommended**: Use [LazyVim](https://lazyvim.org) or `lazy.nvim`.

## ⚡ Accelerated Setup

### 1. Install the Server

```bash
npm install -g dft-language-server
```

This installs the `vasp-lsp` and `dft-lsp` executables.

### 2. Install the Neovim Plugin (Optional)

Using `lazy.nvim`:

```lua
return {
  {
    "TheFermiSea/vasp-language-server",
    ft = { "vasp", "crystal" },
    config = function()
      require("vasp-lsp").setup()
    end,
  },
}
```

### 3. Configure LazyVim

Add this to `lua/plugins/vasp.lua`:

```lua
return {
  {
    "neovim/nvim-lspconfig",
    opts = {
      servers = {
        vasp_ls = {
          cmd = { "vasp-lsp", "--stdio" },
          filetypes = { "vasp", "crystal" },
        },
      },
      setup = {
        vasp_ls = function()
          vim.filetype.add({
            pattern = {
              ['.*INCAR.*'] = 'vasp',
              ['.*POSCAR.*'] = 'vasp',
              ['.*CONTCAR.*'] = 'vasp',
              ['.*POTCAR.*'] = 'vasp',
              ['.*KPOINTS.*'] = 'vasp',
              ['.*%.d12'] = 'crystal',
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

1. **[nvim-cmp](https://github.com/hrsh7th/nvim-cmp)**: Autocompletion for INCAR/CRYSTAL keywords.
2. **[telescope.nvim](https://github.com/nvim-telescope/telescope.nvim)**: `:Telescope lsp_document_symbols` for large files.
3. **[catppuccin](https://github.com/catppuccin/nvim)**: Semantic token highlighting support.

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
    ['.*%.d12'] = 'crystal',
  },
})

-- 2. Setup LSP
local lspconfig = require('lspconfig')
local configs = require('lspconfig.configs')

if not configs.vasp_ls then
  configs.vasp_ls = {
    default_config = {
      cmd = { 'vasp-lsp', '--stdio' },
      filetypes = { 'vasp', 'crystal' },
      root_dir = lspconfig.util.root_pattern('.git', 'INCAR', 'POSCAR'),
      settings = {},
    },
  }
end

lspconfig.vasp_ls.setup{}
```

### Optional: Filetype Overrides

If your files use custom names or extensions, override detection with server settings:

```lua
lspconfig.vasp_ls.setup{
  settings = {
    vasp = {
      fileTypeOverrides = {
        filenames = {
          ["INCAR.relax"] = "incar",
          ["POSCAR.start"] = "poscar",
        },
        extensions = {
          [".vasp"] = "poscar",
          ["kp"] = "kpoints",
        },
      },
    },
  },
}
```
