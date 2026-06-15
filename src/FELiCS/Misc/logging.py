#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
# Standard libraries
from datetime import datetime
import logging
import os
import re
import shutil
import subprocess
import sys

def log_and_raise(logger, log_message, exc_type=RuntimeError, raise_message=None):
    logger.error(log_message)
    raise exc_type(raise_message if raise_message is not None else log_message)

def in_notebook():
    try:
        from IPython import get_ipython
        return get_ipython() is not None
    except ImportError:
        return False

class CustomFormatter(logging.Formatter):
    """
    Custom log formatter with colored output for different log levels.

    Provides colored formatting for log messages in the terminal, making it easier to distinguish between log levels such as DEBUG, INFO,
    WARNING, ERROR, and CRITICAL.

    **Initialize the CustomFormatter object**

    No parameters are required for initialization.

    Attributes
    ----------
    grey : str
        ANSI escape code for grey color.
    yellow : str
        ANSI escape code for yellow color.
    red : str
        ANSI escape code for red color.
    bold_red : str
        ANSI escape code for bold red color.
    reset : str
        ANSI escape code to reset color.
    format : str
        Log message format string.
    FORMATS : dict
        Mapping of log levels to their respective colored format strings.
    """
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    _fmt      = '%(levelname)-8s | %(filename)-22s | %(funcName)-26s (line %(lineno)-4s) : %(message)s'
    _fmt_debug = '%(asctime)s | %(levelname)-8s | %(filename)-22s | %(funcName)-26s (line %(lineno)-4s) : %(message)s'

    def __init__(self, debug_mode=False):
        super().__init__()
        fmt = self._fmt_debug if debug_mode else self._fmt
        self.FORMATS = {
            logging.DEBUG:    self.grey     + fmt + self.reset,
            logging.INFO:     self.grey     + fmt + self.reset,
            logging.WARNING:  self.yellow   + fmt + self.reset,
            logging.ERROR:    self.red      + fmt + self.reset,
            logging.CRITICAL: self.bold_red + fmt + self.reset,
        }
        self._datefmt = '%Y-%m-%d %H:%M:%S' if debug_mode else None

    def format(
        self,
        record,
    ):
        """
        Format the specified log record as text with color based on log level.

        Applies color formatting to the log message depending on the log level.

        Parameters
        ----------
        record : logging.LogRecord
            The log record to be formatted.

        Returns
        -------
        str
            The formatted log message string.
        """
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt=self._datefmt)
        if in_notebook():

            # Strip ANSI codes inside notebooks
            return re.sub(
                r'\x1b\[[0-9;]*m',
                '',
                formatter.format(record),
            )
        else:
            return formatter.format(record)

class Logger:
    """
    Logger class for flexible and colored logging to file and console.

    Provides a singleton logger with support for colored console output, file logging, dynamic log file location changes, and debug/test modes.
    **Initialize the Logger object**

    Parameters
    ----------
    debug_mode : bool, optional
        If True, enables debug logging (default is False).
    test_mode : bool, optional
        If True, enables test mode logging (default is False).
    profiler_mode : bool, optional
        If True, enables test mode logging (default is False).
    logger_name : str, optional
        Name of the logger and log file prefix (default is "log").

    Attributes
    ----------
    debug_mode : bool
        Indicates if debug mode is enabled.
    test_mode : bool
        Indicates if test mode is enabled.
    profiler_mode : bool
        Indicates if profilers are enabled.
    logger_name : str
        Name of the logger.
    _logger : logging.Logger
        The underlying Python logger instance.
    """

    _instance = None
    
    def __init__(
        self,
        debug_mode=False,
        test_mode=False,
        profiler_mode=False,
        logger_name="log"
    ):
        """
        Initialize the Logger instance.

        Parameters
        ----------
        debug_mode : bool, optional
            If True, enables debug logging (default is False).
        test_mode : bool, optional
            If True, enables test mode logging (default is False).
        profiler_mode : bool, optional
            If True, enables profilers (default is False).
        logger_name : str, optional
            Name of the logger and log file prefix (default is "log").
        """
        self.debug_mode = debug_mode
        self.test_mode = test_mode
        self.profiler_mode = profiler_mode
        self.logger_name = logger_name
        self._logger = None
        self._setup_logger()

        Logger._instance = self._logger

    def _setup_logger(
        self
    ):
        """
        Set up the logger with appropriate handlers and formatters.

        Configures the logger to output to both the console (with colored output)
        and a log file. Creates the logs directory if it does not exist.
        """
        now = datetime.now().strftime("%d.%m.%Y-%H.%M.%S.%f")
        
        logfilename = f"logs{os.sep}{self.logger_name}_{now}.log"

        #logfilename_errors = f"logs{os.sep}{self.logger_name}_{now}.ERRORS.log"

        self._logger = logging.getLogger(self.logger_name)
        self._logger.handlers.clear()
        self._logger.propagate = False
        self._logger._felics_profiler_mode = self.profiler_mode
        
        self._logger.setLevel(logging.DEBUG if (self.debug_mode or self.test_mode) else logging.INFO)

        if not os.path.exists("logs"):
            os.makedirs("logs")

        ch = logging.StreamHandler()
        fh = logging.FileHandler(
            logfilename,
            encoding='utf-8',
        )

        #fh_errors = logging.FileHandler(logfilename_errors, encoding='utf-8')

        ch.setLevel(logging.ERROR if self.test_mode else 
                   logging.DEBUG if self.debug_mode else logging.INFO)
        fh.setLevel(logging.DEBUG if (self.debug_mode or self.test_mode) else logging.INFO)

        #fh_errors.setLevel(logging.WARNING)

        logging.addLevelName(
            logging.ERROR,
            'Error',
        )
        logging.addLevelName(
            logging.WARNING,
            'Warning',
        )
        logging.addLevelName(
            logging.INFO,
            'Info',
        )
        logging.addLevelName(
            logging.DEBUG,
            'Debug',
        )
        _log_fmt = (
            '%(asctime)s | %(levelname)-8s | %(filename)-22s | %(funcName)-26s (line %(lineno)-4s) : %(message)s'
            if self.debug_mode else
            '%(levelname)-8s | %(filename)-22s | %(funcName)-26s (line %(lineno)-4s) : %(message)s'
        )
        formatter_log = logging.Formatter(
            _log_fmt,
            datefmt='%Y-%m-%d %H:%M:%S' if self.debug_mode else None,
        )

        # use custom formatter to get colored output for the command line
        formatter_cmd = CustomFormatter(debug_mode=self.debug_mode)
        
        ch.setFormatter(formatter_cmd)
        fh.setFormatter(formatter_log)

        #fh_errors.setFormatter(formatter_log)
        for handler in [ch, fh]: #, fh_errors]:
            self._logger.addHandler(handler)

        self._logger.debug("Logger initialized successfully")
        
        if self.profiler_mode:
            self._start_mem_tracker()
            self._start_cpu_tracker()

    def _start_cpu_tracker(self):
        """
        Launch the background CPU-usage tracker alongside the log file.

        Resolves the path to tools/cpu_tracker.py relative to this file,
        then spawns it as a detached subprocess that samples CPU utilization and writes
        a CSV and PNG into the same directory as the log file.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        log_file = next(
            (h.baseFilename for h in self._logger.handlers
             if isinstance(h, logging.FileHandler)),
            None,
        )
        if log_file is None:
            return

        tracker = os.path.normpath(
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', '..', '..', 'tools', 'cpu_tracker.py')
        )
        if not os.path.isfile(tracker):
            self._logger.warning(f"CPU tracker script not found: {tracker}")
            return

        subprocess.Popen(
            [sys.executable, tracker,
             '--pid',      str(os.getpid()),
             '--log-file', log_file],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        self._logger.info(
            f"CPU tracker started for PID {os.getpid()}; "
            f"output in log file directory"
        )

    def _start_mem_tracker(self):
        """
        Launch the background memory-usage tracker alongside the log file.

        Resolves the path to tools/memory_tracker.py relative to this file,
        then spawns it as a detached subprocess that samples memory and writes
        a CSV and PNG into the same directory as the log file.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        log_file = next(
            (h.baseFilename for h in self._logger.handlers
             if isinstance(h, logging.FileHandler)),
            None,
        )
        if log_file is None:
            return

        tracker = os.path.normpath(
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', '..', '..', 'tools', 'memory_tracker.py')
        )
        if not os.path.isfile(tracker):
            self._logger.warning(f"Memory tracker script not found: {tracker}")
            return

        subprocess.Popen(
            [sys.executable, tracker,
             '--pid',      str(os.getpid()),
             '--log-file', log_file],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        self._logger.info(
            f"Memory tracker started for PID {os.getpid()}; "
            f"output in log file directory"
        )

    @property
    def logger(
        self
    ):
        """
        Get the underlying logger instance.

        Returns
        -------
        logging.Logger
            The configured logger instance.

        Raises
        ------
        RuntimeError
            If the logger has not been properly initialized.
        """
        if self._logger is None:
            raise RuntimeError("Logger has not been properly initialized")
        return self._logger

    @classmethod
    def get_logger(
        cls,
        name=None
    ):
        """
        Retrieve the singleton logger instance.

        If the logger has not been initialized, creates a new instance with the specified name.

        Parameters
        ----------
        name : str, optional
            Name for the logger if it needs to be created.

        Returns
        -------
        logging.Logger
            The singleton logger instance.

        Raises
        ------
        RuntimeError
            If the logger is not initialized and no name is provided.
        """
        if cls._instance is None:
            if name is None:
                raise RuntimeError("Logger not initialized. Create a Logger instance first.")
            return cls(logger_name=name).logger
        return cls._instance
    
    @classmethod
    def change_log_location(
        cls,
        new_log_path
    ):
        """
        Move log files to a new directory and update file handlers.

        Moves all current log files to the specified new directory and updates the logger's file handlers to write to the new location.

        Renames memory-profiler companion files (CSV, PNG) if they exist, and writes a redirect file so the background tracker switches its output directory as well.

        Parameters
        ----------
        new_log_path : str
            The new directory path for log files.

        Raises
        ------
        RuntimeError
            If the logger is not initialized.
        """
        if cls._instance is None:
            raise RuntimeError("Logger not initialized. Create a Logger instance first.")

        instance = cls.get_logger()

        os.makedirs(
            os.path.dirname(new_log_path),
            exist_ok=True,
        )

        file_handlers = [h for h in instance.handlers if isinstance(
            h,
            logging.FileHandler,
        )]

        for handler in file_handlers:
            log_level = handler.level
            log_formatter = handler.formatter

            instance.removeHandler(handler)
            handler.close()

            if not os.path.exists(new_log_path):
                os.makedirs(new_log_path)

            old_log_path = handler.baseFilename
            if os.path.exists(old_log_path):
                shutil.move(
                    old_log_path,
                    new_log_path,
                )

            # Move tracker files (CSV, PNG) 
            # if they exist, write a redirect file so the background tracker switches its output directory.
            old_dir  = os.path.dirname(old_log_path)
            log_stem = os.path.splitext(os.path.basename(old_log_path))[0]
            profiler_active = bool(getattr(instance, "_felics_profiler_mode", False))
            for _tracker in ["memory", "cpu"]:
                _has_tracker_files = False
                for _ext in ["csv", "png"]:
                    _src = os.path.join(old_dir, f"{log_stem}_{_tracker}.{_ext}")
                    if os.path.exists(_src):
                        shutil.move(_src, new_log_path)
                        _has_tracker_files = True

                if _has_tracker_files or profiler_active:
                    _redirect = os.path.join(old_dir, f"{log_stem}_{_tracker}.redirect")
                    try:
                        with open(_redirect, "w") as _rf:
                            _rf.write(new_log_path)
                    except OSError:
                        pass
            
            new_log_file_path = os.path.join(
                new_log_path,
                os.path.basename(old_log_path),
            )

            new_file_handler = logging.FileHandler(
                new_log_file_path,
                encoding='utf-8',
            )
            new_file_handler.setLevel(log_level)
            new_file_handler.setFormatter(log_formatter)

            instance.addHandler(new_file_handler)

        instance.info(f"Log files moved to: {new_log_path}")

    def close_logger(
        self,
        logger
    ):
        """
        Remove all handlers from the given logger.

        Clears all handlers from the specified logger instance.

        Parameters
        ----------
        logger : logging.Logger
            The logger instance to close.
        """
        if logger.hasHandlers():
            logger.handlers.clear()
