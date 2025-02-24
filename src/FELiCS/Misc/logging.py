import os
import logging
import shutil
from datetime import datetime

class Logger:
    # Class variable to store the logger instance
    _instance = None
    
    def __init__(self, debug_mode=False, logger_name="log"):
        self.debug_mode     = debug_mode
        self.logger_name    = logger_name
        self._logger        = None
        self._setup_logger()

        # Store logger instance in class variable
        Logger._instance = self._logger

    def _setup_logger(self):
        # Get time stamp
        now                 = datetime.now().strftime("%d.%m.%Y-%H.%M.%S")
        
        # Define log files
        logfilename         = f"logs{os.sep}{self.logger_name}_{now}.log"
        logfilename_errors  = f"logs{os.sep}{self.logger_name}_{now}.ERRORS.log"

        # Create logger
        self._logger        = logging.getLogger(self.logger_name)
        
        # Clear any existing handlers and prevent propagation to root logger
        self._logger.handlers.clear()
        self._logger.propagate = False
        
        self._logger.setLevel(logging.DEBUG if self.debug_mode else logging.INFO)

        # Create logs directory
        if not os.path.exists("logs"):
            os.makedirs("logs")

        # Setup handlers
        ch          = logging.StreamHandler()
        fh          = logging.FileHandler(logfilename, encoding='utf-8')
        fh_errors   = logging.FileHandler(logfilename_errors, encoding='utf-8')

        # Set levels
        ch.setLevel(logging.DEBUG if self.debug_mode else logging.INFO)
        fh.setLevel(logging.DEBUG if self.debug_mode else logging.INFO)
        fh_errors.setLevel(logging.WARNING)

        # Setup formatter
        logging.addLevelName(logging.ERROR,     'Error')
        logging.addLevelName(logging.WARNING,   'Warning')
        logging.addLevelName(logging.INFO,      'Info')
        logging.addLevelName(logging.DEBUG,     'Debug')
        formatter = logging.Formatter(
            '%(levelname)-8s: %(filename)-16s | %(funcName)-24s (line %(lineno)-4s) : %(message)s'
        )

        # Apply formatter
        for handler in [ch, fh, fh_errors]:
            handler.setFormatter(formatter)
            self._logger.addHandler(handler)

        self._logger.debug("Logger initialized successfully")

    @property
    def logger(self):
        """Get the logger instance."""
        if self._logger is None:
            raise RuntimeError("Logger has not been properly initialized")
        return self._logger

    @classmethod
    def get_logger(cls, name=None):
        """Get the global logger instance.
        
        Args:
            name (str, optional): Logger name if creating new logger. Defaults to None.
        
        Returns:
            logging.Logger: Logger instance
        """
        if cls._instance is None:
            if name is None:
                raise RuntimeError("Logger not initialized. Create a Logger instance first.")
            return cls(logger_name=name).logger
        return cls._instance
    
    @classmethod
    def change_log_location(cls, new_log_path):
        """Change log file location for the global logger instance.
        
        Args:
            new_log_path (str): New path for log files
        """
        if cls._instance is None:
            raise RuntimeError("Logger not initialized. Create a Logger instance first.")
        
        instance = cls.get_logger()
        
        # Ensure the new directory exists
        os.makedirs(os.path.dirname(new_log_path), exist_ok=True)
        
        # Get the current file handlers
        file_handlers = [h for h in instance.handlers if isinstance(h, logging.FileHandler)]
        
        for handler in file_handlers:
            # Store existing log level and formatter
            log_level       = handler.level
            log_formatter   = handler.formatter

            # Close and remove the existing file handler
            instance.removeHandler(handler)
            handler.close()
            
            # Create new log folder if needed
            if not os.path.exists(new_log_path):
                os.makedirs(new_log_path)

            # Move the old log file to the new location (if it exists)
            old_log_path = handler.baseFilename
            if os.path.exists(old_log_path):
                shutil.move(old_log_path, new_log_path)
            new_log_file_path = os.path.join(new_log_path, os.path.basename(old_log_path))
            
            # Create a new file handler with the same settings
            new_file_handler = logging.FileHandler(new_log_file_path, encoding='utf-8')
            new_file_handler.setLevel(log_level)            # Preserve log level
            new_file_handler.setFormatter(log_formatter)    # Preserve formatter

            # Attach the new handler to the logger
            instance.addHandler(new_file_handler)

        # Log the change
        instance.info(f"Log files moved to: {new_log_path}")

    def close_logger(self, logger):
        """
        Parameters
        ----------
        logger : _type_
            _description_
        """
        if logger.hasHandlers():
            logger.handlers.clear()

