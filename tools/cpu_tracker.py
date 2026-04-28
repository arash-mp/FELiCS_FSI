#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \  |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/  |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
# Background CPU utilization tracker for FELiCS runs.
#
# Reads CPU time metrics from /proc/<pid>/stat (Linux only),
# writes a CSV, and regenerates a PNG plot after each sample.
# Only standard-library modules plus matplotlib are used.
#
# Usage (called automatically by logger in debug mode):
#   python cpu_tracker.py --pid <PID> --log-file <path/to/felics_<timestamp>.log>
#                             [--interval <seconds>]

import argparse
import csv
import os
import time
from datetime import datetime, timedelta

import matplotlib
matplotlib.use('Agg')           # non-interactive backend – no display needed
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Suffix used to communicate a log-directory change from Logger.change_log_location
_REDIRECT_SUFFIX = "_cpu.redirect"


# ---------------------------------------------------------------------------
# /proc helpers (Linux)
# ---------------------------------------------------------------------------

def _read_cpu_ticks(pid):
    """
    Return the total CPU ticks (user + system time) for `pid` by parsing /proc/<pid>/stat.

    Parameters
    ----------
    pid : int
        Process ID to read CPU status for.

    Returns
    -------
    int or None
        Total CPU ticks so far, or None if the process no longer exists.
    """
    stat_path = f"/proc/{pid}/stat"
    try:
        with open(stat_path, 'r') as fh:
            data = fh.read()
            
        # Structure of /proc/[pid]/stat has the command name in parentheses.
        # This name can contain spaces, so partition by the last closing parenthesis.
        pid_s, rest = data.split('(', 1)
        comm, rest = rest.rsplit(')', 1)
        fields = rest.split()
        
        # fields[0] corresponds to state, which is the 3rd field in the stat file overall
        # utime is 14th field -> index 11 in `fields`
        # stime is 15th field -> index 12 in `fields`
        utime = int(fields[11])
        stime = int(fields[12])
        
        return utime + stime
    except (OSError, ValueError, IndexError):
        return None


def _process_exists(pid):
    """
    Check if a process with the given PID exists.

    Parameters
    ----------
    pid : int
        Process ID to check.

    Returns
    -------
    bool
        True if the process exists in /proc, False otherwise.
    """
    return os.path.exists(f"/proc/{pid}")


# ---------------------------------------------------------------------------
# Plot helper
# ---------------------------------------------------------------------------

def _update_plot(timestamps, cpu_percentages, png_path, elapsed_times):
    """
    Generate and save a diagnostic plot of CPU utilization over time.

    Parameters
    ----------
    timestamps : list of datetime.datetime
        List of timestamps corresponding to each CPU sample.
    cpu_percentages : list of float
        List of CPU utilization percentages.
    png_path : str
        File path where the generated plot image will be saved.
    elapsed_times : list of float
        List of elapsed times from the start of the tracking.

    Returns
    -------
    None
    """
    max_cpu = max(cpu_percentages) if cpu_percentages else 0.0
    total_runtime = elapsed_times[-1] if elapsed_times else 0.0

    # Format total runtime efficiently
    runtime_str = str(timedelta(seconds=int(total_runtime))) if total_runtime >= 60 else f"{total_runtime:.1f}s"

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(timestamps, cpu_percentages,
            label=f'CPU (max: {max_cpu:.1f}%)', color='tab:red')
    ax.set_xlabel('Time')
    ax.set_ylabel('CPU Utilization (%)')
    ax.set_title(f'FELiCS CPU Usage based on CPU clock ticks per second (Total Runtime: {runtime_str})')
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    fig.autofmt_xdate()
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(png_path, dpi=100)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """
    Main entry point for the CPU tracker script.

    Parses command-line arguments to monitor the CPU usage of a specified PID,
    writes the tracking data to a CSV file, and iteratively updates a plot image.
    Outputs are placed alongside the provided log file path and track redirects
    if the log file location changes during execution.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(
        description='Track CPU utilization of a FELiCS process and '
                    'store results alongside the log file.'
    )
    parser.add_argument('--pid', type=int, required=True,
                        help='PID of the process to monitor')
    parser.add_argument('--log-file', type=str, required=True,
                        help='Path to the FELiCS log file '
                             '(CSV and PNG are placed in the same directory)')
    parser.add_argument('--interval', type=float, default=1.0,
                        help='Sampling interval in seconds (default: 1.0)')
    args = parser.parse_args()

    # Determine CPU clock ticks per second (usually 100 on Linux)
    try:
        clk_tck = os.sysconf("SC_CLK_TCK")
    except (ValueError, AttributeError):
        clk_tck = 100.0

    # Derive output paths from the log file path
    log_dir  = os.path.dirname(os.path.abspath(args.log_file))
    log_stem = os.path.splitext(os.path.basename(args.log_file))[0]
    csv_path = os.path.join(log_dir, f"{log_stem}_cpu.csv")
    png_path = os.path.join(log_dir, f"{log_stem}_cpu.png")

    os.makedirs(log_dir, exist_ok=True)

    # Write CSV header
    with open(csv_path, 'w', newline='') as fh:
        csv.writer(fh).writerow(['timestamp', 'elapsed_s', 'cpu_percent'])

    elapsed_times   = []
    timestamps      = []
    cpu_percentages = []
    
    start_time = time.time()
    
    # Initialize previous state for CPU delta calculation
    prev_time = start_time
    prev_ticks = _read_cpu_ticks(args.pid)
    
    if prev_ticks is None:
        # Cannot read CPU info, perhaps process died immediately
        return

    while _process_exists(args.pid):
        time.sleep(args.interval)
        
        # Check whether Logger.change_log_location has moved the log file.
        redirect_file = os.path.join(log_dir, f"{log_stem}{_REDIRECT_SUFFIX}")
        if os.path.exists(redirect_file):
            try:
                with open(redirect_file) as rf:
                    new_dir = rf.read().strip()
                os.remove(redirect_file)
                log_dir  = new_dir
                csv_path = os.path.join(log_dir, f"{log_stem}_cpu.csv")
                png_path = os.path.join(log_dir, f"{log_stem}_cpu.png")
                
                if not os.path.exists(csv_path):
                    os.makedirs(log_dir, exist_ok=True)
                    with open(csv_path, 'w', newline='') as fh:
                        csv.writer(fh).writerow(
                            ['timestamp', 'elapsed_s', 'cpu_percent']
                        )
            except OSError:
                pass

        current_ticks = _read_cpu_ticks(args.pid)
        if current_ticks is None:
            break

        current_time = time.time()
        
        delta_ticks = current_ticks - prev_ticks
        delta_time = current_time - prev_time
        
        cpu_percent = 0.0
        if delta_time > 0:
            # (ticks / ticks_per_second) / total_elapsed_seconds
            cpu_percent = 100.0 * (delta_ticks / clk_tck) / delta_time

        now     = datetime.now()
        elapsed = current_time - start_time
        ts      = now.strftime("%H:%M:%S.%f")

        elapsed_times.append(elapsed)
        timestamps.append(now)
        cpu_percentages.append(cpu_percent)

        # Append row to CSV
        with open(csv_path, 'a', newline='') as fh:
            csv.writer(fh).writerow(
                [ts, f"{elapsed:.3f}", f"{cpu_percent:.1f}"]
            )

        # Regenerate plot
        _update_plot(timestamps, cpu_percentages, png_path, elapsed_times)

        # Slide window forward
        prev_ticks = current_ticks
        prev_time = current_time


if __name__ == '__main__':
    main()