-- vasp-lsp.lua
-- Unified configuration for the VASP Language Server in Neovim.

local M = {}

--- Registers the VASP filetypes.
function M.register_filetypes()
  vim.filetype.add({
    pattern = {
      ['.*INCAR.*'] = 'vasp',
      ['.*POSCAR.*'] = 'vasp',
      ['.*CONTCAR.*'] = 'vasp',
      ['.*POTCAR.*'] = 'vasp',
      ['.*KPOINTS.*'] = 'vasp',
    },
  })
end

--- Returns the default configuration for lspconfig.
--- @param cmd? table The command to run the server. Defaults to {'vasp-lsp', '--stdio'} if found in path, else node /path/to/server.
function M.get_default_config(cmd)
  local lspconfig = require('lspconfig')
  
  return {
    cmd = cmd or { 'vasp-lsp', '--stdio' },
    filetypes = { 'vasp' },
    root_dir = lspconfig.util.root_pattern('.git', 'INCAR', 'POSCAR'),
    settings = {},
    -- Default semantic tokens legend for VASP
    on_attach = function(client, _)
      if client.server_capabilities.semanticTokensProvider then
        client.server_capabilities.semanticTokensProvider.legend = {
          tokenTypes = { "property", "keyword", "number", "comment", "string" },
          tokenModifiers = { "declaration" }
        }
      end
    end
  }
end

--- Easy setup function for lspconfig.
--- @param opts? table Custom options to pass to setup.
function M.setup(opts)
  M.register_filetypes()
  
  local lspconfig = require('lspconfig')
  local configs = require('lspconfig.configs')
  
  if not configs.vasp_ls then
    configs.vasp_ls = {
      default_config = M.get_default_config(opts and opts.cmd)
    }
  end
  
  lspconfig.vasp_ls.setup(opts or {})
end

return M
