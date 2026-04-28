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
# Background memory tracker for FELiCS runs.
#
# Reads resident/virtual memory from /proc/<pid>/status (Linux only),
# writes a CSV, and regenerates a PNG plot after each sample.
# Only standard-library modules plus matplotlib (already in the felics env) are used.
#
# Usage (called automatically by main.py in debug mode):
#   python track_mem_usage.py --pid <PID> --log-file <path/to/felics_<timestamp>.log>
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
_REDIRECT_SUFFIX = "_memory.redirect"


# ---------------------------------------------------------------------------
# /proc helpers (Linux)
# ---------------------------------------------------------------------------

def _read_proc_mem_kb(pid):
    """
    Return (rss_kb, vsz_kb) for `pid` by parsing /proc/<pid>/status.

    Parameters
    ----------
    pid : int
        Process ID to read memory status for.

    Returns
    -------
    tuple of (int or None, int or None)
        A tuple containing Resident Set Size (RSS) and Virtual Memory Size (VSZ) in kilobytes.
        Returns (None, None) if the process no longer exists or the file cannot be read.
    """
    status_path = f"/proc/{pid}/status"
    rss_kb = vsz_kb = None
    try:
        with open(status_path, 'r') as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    rss_kb = int(line.split()[1])
                elif line.startswith("VmSize:"):
                    vsz_kb = int(line.split()[1])
                if rss_kb is not None and vsz_kb is not None:
                    break
    except OSError:
        return None, None
    return rss_kb, vsz_kb


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

def _update_plot(timestamps, rss_mbs, vsz_mbs, png_path):
    """
    Generate and save a diagnostic plot of memory usage over time.

    Parameters
    ----------
    timestamps : list of datetime.datetime
        List of timestamps corresponding to each memory sample.
    rss_mbs : list of float
        List of Resident Set Size memory values in megabytes.
    vsz_mbs : list of float
        List of Virtual Memory Size values in megabytes.
    png_path : str
        File path where the generated plot image will be saved.

    Returns
    -------
    None
    """
    rss_gbs = [v / 1024.0 for v in rss_mbs]
    vsz_gbs = [v / 1024.0 for v in vsz_mbs]
    max_rss = max(rss_gbs)
    max_vsz = max(vsz_gbs)

    # Format total runtime efficiently
    total_runtime = (timestamps[-1] - timestamps[0]).total_seconds() if timestamps else 0.0
    runtime_str = str(timedelta(seconds=int(total_runtime))) if total_runtime >= 60 else f"{total_runtime:.1f}s"

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(timestamps, rss_gbs,
            label=f'RSS (max: {max_rss:.2f} GB)', color='tab:blue') 
    ax.plot(timestamps, vsz_gbs,
            label=f'Virtual (max: {max_vsz:.2f} GB)',
            color='tab:orange', linestyle='--')
    ax.set_xlabel('Time')
    ax.set_ylabel('Memory (GB)')
    ax.set_title(f'FELiCS Memory Usage (Total Runtime: {runtime_str})')
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
    Main entry point for the memory tracker script.

    Parses command-line arguments to monitor the memory usage of a specified PID,
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
        description='Track resident/virtual memory of a FELiCS process and '
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

    # Derive output paths from the log file path
    log_dir  = os.path.dirname(os.path.abspath(args.log_file))
    log_stem = os.path.splitext(os.path.basename(args.log_file))[0]
    csv_path = os.path.join(log_dir, f"{log_stem}_memory.csv")
    png_path = os.path.join(log_dir, f"{log_stem}_memory.png")

    os.makedirs(log_dir, exist_ok=True)

    # Write CSV header
    with open(csv_path, 'w', newline='') as fh:
        csv.writer(fh).writerow(['timestamp', 'elapsed_s', 'rss_mb', 'vsz_mb'])

    elapsed_times = []
    timestamps    = []
    rss_mbs       = []
    vsz_mbs       = []
    start_time    = time.time()

    while _process_exists(args.pid):
        # Check whether Logger.change_log_location has moved the log file.
        # If so, a redirect file is left in the old directory containing the
        # new directory path.  Follow it so CSV and PNG end up alongside the log.
        redirect_file = os.path.join(log_dir, f"{log_stem}{_REDIRECT_SUFFIX}")
        if os.path.exists(redirect_file):
            try:
                with open(redirect_file) as rf:
                    new_dir = rf.read().strip()
                os.remove(redirect_file)
                log_dir  = new_dir
                csv_path = os.path.join(log_dir, f"{log_stem}_memory.csv")
                png_path = os.path.join(log_dir, f"{log_stem}_memory.png")
                # Companion files were moved by change_log_location; write a
                # fresh CSV header only if the moved file is not already there.
                if not os.path.exists(csv_path):
                    os.makedirs(log_dir, exist_ok=True)
                    with open(csv_path, 'w', newline='') as fh:
                        csv.writer(fh).writerow(
                            ['timestamp', 'elapsed_s', 'rss_mb', 'vsz_mb']
                        )
            except OSError:
                pass

        rss_kb, vsz_kb = _read_proc_mem_kb(args.pid)
        if rss_kb is None:
            break

        now     = datetime.now()
        elapsed = time.time() - start_time
        ts      = now.strftime("%H:%M:%S.%f")
        rss_mb  = rss_kb / 1024.0
        vsz_mb  = vsz_kb / 1024.0

        elapsed_times.append(elapsed)
        timestamps.append(now)
        rss_mbs.append(rss_mb)
        vsz_mbs.append(vsz_mb)

        # Append row to CSV
        with open(csv_path, 'a', newline='') as fh:
            csv.writer(fh).writerow(
                [ts, f"{elapsed:.3f}", f"{rss_mb:.2f}", f"{vsz_mb:.2f}"]
            )

        # Regenerate plot
        _update_plot(timestamps, rss_mbs, vsz_mbs, png_path)

        time.sleep(args.interval)


if __name__ == '__main__':
    main()
