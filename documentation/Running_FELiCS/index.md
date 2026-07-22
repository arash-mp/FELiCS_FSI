# Running-FELiCS

The typical *monolithic* FELiCS workflow consists of:
1. Preparing the required [input files](#input-files)
2. Running the [`FELiCS` command](#felics-commands)
3. Post-processing the generated [output files](#output-files)

(input-files)=
## Input files

The main input files for FELiCS are:
- a configuration file (`.json`)
- a mesh file (`.msh`)
- a boundary conditions file (`.json`)
- a mean flow file (`.fel` or `MeanFlow.h5`)

Additional input files may be required depending on the simulation setup.

The different input files are described in the following sections:

| File | Description |
| - | - |
| [Configuration file](input_files/param_files.md) | Main simulation settings and solver parameters |
| [Mesh file](input_files/gmsh_files.md) | Computational mesh in Gmsh format |
| [Mean flow file](input_files/mean_files.md) | Mean flow solution used by FELiCS |
| [Boundary file](input_files/bc_files.md) | Boundary condition definitions |
| [Mixture file](input_files/mix_files.md) | Optional mixture properties for reacting flows |

(felics-commands)=
## FELiCS Commands

After completing the steps in the [installation guide](../installation_guide.md), the `FELiCS` program can be called directly from the command line, for instance with:

```bash
FELiCS -f <your_config_file.json>
```

For an overview of the arguments for the `FELiCS` command:

```bash
FELiCS -h
```

| Argument | Description |
| - | - |
| `-h` | Show the help message and exit |
| `-f <path>`, `--file <path>` | (**Mandatory**) Specifies the path to a FELiCS [config file](input_files/param_files.md) (in `json` format). |
| `-d`, `--debug` | (*Optional*) Enable *debug mode* for extended verbose output |
| `-t`, `--test` | (*Optional*) Enable *test mode* with no verbose output |

Example usage:

```bash
FELiCS -f config.json
FELiCS -f config.json -d
FELiCS -f config.json -t
```

(output-files)=
## Output files

After running an analysis, FELiCS generates several output files containing:
- computed solutions (`.h5`),
- interpolated mean-flow fields (`.h5`),
- computational meshes (`.h5`),
- solution spectra (`.csv`),
- and execution logs (`.log`).

For visualization purposes, FELiCS also automatically generates associated `.xmf` metadata files compatible with ParaView.

The generated files depend on the selected analysis type and export options.

The different output files are described in the following sections:

| File | Description |
| - | - |
| [Solution files](output_files/solution_files.md) | Modal, resolvent, or input/output solutions exported by FELiCS |
| [FELiCS Mean-flow files](output_files/FELiCSmeanflow_files.md) | Mean-flow fields exported onto the FELiCS mesh |
| [FELiCS Mesh files](output_files/FELiCSmesh_files.md) | Exported FEM meshes used by FELiCS |
| [Spectrum files](output_files/spectrum_files.md) | Eigenvalue or gain spectra exported as CSV files |
| [Log files](output_files/log_files.md) | Solver logs and execution information |

## Summary

The diagram below summarizes the main inputs and outputs of FELiCS. It also provides links to the corresponding documentation pages and tutorials.

```{mermaid}
flowchart LR
    subgraph Input["Input Data"]
        A[Mean or Base<br>Flow Data<br> <a href="../Tutorials/cylinder_wake.html">see Tutorial 1</a>]
        click A "../Tutorials/cylinder_wake.html" _blank
        style Input fill:#f1faee,stroke:#457b9d,stroke-width:2px
        style A fill:#cdb4db,stroke:#333,stroke-width:1px
    end

    subgraph PrepFiles["FELiCS input files"]
        B1[settings.json<br><a href="input_files/param_files.html">see details</a>]
        B2[.fel file<br><a href="input_files/mean_files.html">see details</a>]
        B3[.msh file<br><a href="../Tutorials/cylinder_wake.html">see Tutorial 1</a>]
        B4[boundaries.json<br><a href="input_files/bc_files.html">see details</a><br><a href="../Tutorials/modal_analysis.html">see Tutorial 2</a>]
        B5[mixture.json<br><a href="input_files/mix_files.html">see details</a>]

        style PrepFiles fill:#f1faee,stroke:#457b9d,stroke-width:2px
        style B1 fill:#ffafcc,stroke:#333,stroke-width:1px
        style B2 fill:#ffafcc,stroke:#333,stroke-width:1px
        style B3 fill:#ffafcc,stroke:#333,stroke-width:1px
        style B4 fill:#ffafcc,stroke:#333,stroke-width:1px
        style B5 fill:#ffafcc,stroke:#333,stroke-width:1px


    end

    C[FELiCS<br><a href="#felics-commands">see details</a>]
    style C fill:#e76f51

    subgraph Output["Output"]
        D1[Modal Analysis<br><a href="../Tutorials/modal_analysis.html">see Tutorial 2</a>]
        D2[Resolvent Analysis<br><a href="../Tutorials/Resolvent_Analysis.html">see Tutorial 3</a>]
        D3[Input-Output<br>Analysis<br><a href="../Tutorials/input_output_analysis.html">see Tutorial 4</a>]
        style Output fill:#f1faee,stroke:#457b9d,stroke-width:2px
        style D1 fill:#a8dadc,stroke:#333,stroke-width:1px
        style D2 fill:#a8dadc,stroke:#333,stroke-width:1px
        style D3 fill:#a8dadc,stroke:#333,stroke-width:1px


    end

    A ==> B2
    B1 ==> C


    B2 & B3 & B4 & B5 -. path .-> B1

    C ==> D1
    C ==> D2
    C ==> D3

    linkStyle default stroke:#8e8d8f,stroke-width:2px


```

```{toctree}
:maxdepth: 2

input_files/index
output_files/index
```
