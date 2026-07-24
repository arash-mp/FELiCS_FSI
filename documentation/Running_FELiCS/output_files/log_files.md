# Log files

FELiCS automatically exports log files containing execution and solver information.

The log files are stored in a `log/` directory inside the FELiCS output directory defined by [Export.ExportFolder](../input_files/parameters/Export/ExportFolder.md).

Typical file names are:

```text
log/
├── felics_29.11.2025-10.44.56.log
└── felics_29.11.2025-10.45.01.log
```

The verbosity and contents of the log files depend on the [command-line arguments](../index.md#felics-commands) used when running FELiCS.

If the profiler is active (using the '-p' flag) there are additional csv and png files showing the CPU usage and the Memory usage.

```{caution}
The CPU and memory profiler only work on Linux.
```

```text
log/
├── felics_29.11.2025-10.44.56.log
├── felics_29.11.2025-10.45.01.log
├── felics_29.11.2025-10.44.56_cpu.csv
├── felics_29.11.2025-10.45.01_cpu.png
├── felics_29.11.2025-10.45.01_memory.csv
└── felics_29.11.2025-10.45.01_memory.png
```

