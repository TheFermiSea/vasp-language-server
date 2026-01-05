-- vasp-lsp plugin entrypoint for Neovim.

local function setup(opts)
  require('vasp-lsp').setup(opts or {})
end

vim.api.nvim_create_user_command('VaspLspSetup', function()
  setup({})
end, {})

if vim.g.vasp_lsp_auto_setup then
  setup(vim.g.vasp_lsp_config)
end
