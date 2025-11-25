# FELiCS Commands

The `FELiCS` program can be called directly from the command line after completing the steps in the installation guide.

## Overview of Arguments

To see an overview of the arguments for the `FELiCS` command, use:

```bash
FELiCS -h
```

This command will display the following options:

| Argument | Description |
| - | - |
| `-h` | Show the help message and exit |
| `-f <path>` or `--file <path>` | (Mandatory) Specifies the path to a FELiCS config file (in `json` format). For more details, refer to [FELiCS Settings](FELiCS_settings.md). |
| `-d` or `--debug` | (Optional) Activate debug mode for extended verbose output |
| `-t` or `--test` | (Optional) Activate test mode which produces no verbose output |

## Example Usages

Here are some examples of how to use the `FELiCS` command:

| Command | Description |
| - | - |
| `FELiCS -h` | Shows help message |
| `FELiCS -f config.json` | Start with config file |
| `FELiCS -f config.json -d` | Start in debug mode |
| `FELiCS -f config.json -t` | Start in test mode |
