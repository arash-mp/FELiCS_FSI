# FELiCS commands

The `FELiCS` program can be called directly from the command line after completion of the steps in the installation guide.

For an overview of the arguments to the `FELiCS` command, use:

```
FELiCS -h
```

Which will provide the following:
* `-h`: show the help message and exit
* `-f <path>` or `--file <path>`: (mandatory) specifies the path to a FELiCS config file (in `json` format) detailed in [this section](FELiCS_settings.md). 
* `-d` or `--debug`: (*optional*) activate debug mode for extended verbose output
* `-t` or `--test`: (*optional*) activate test mode which produces no verbose output.

Below are some example usages:

```
FELiCS -h                   shows help message
FELiCS -f config.json       start with config file
FELiCS -f config.json -d    for debug mode
FELiCS -f config.json -t    for test mode
```
