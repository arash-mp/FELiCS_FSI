import os
import logging
import shutil
from datetime import datetime



class CustomFormatter(logging.Formatter):
    # This class is created to get colored output for the warnings and errors

    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = '%(levelname)-8s | %(filename)-22s | %(funcName)-26s (line %(lineno)-4s) : %(message)s'


    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


class Logger:
    _instance = None
    
    def __init__(self, debug_mode=False, test_mode=False, logger_name="log"):
        self.debug_mode = debug_mode
        self.test_mode = test_mode
        self.logger_name = logger_name
        self._logger = None
        self._setup_logger()

        Logger._instance = self._logger

    def _setup_logger(self):
        now = datetime.now().strftime("%d.%m.%Y-%H.%M.%S")
        
        logfilename = f"logs{os.sep}{self.logger_name}_{now}.log"
        logfilename_errors = f"logs{os.sep}{self.logger_name}_{now}.ERRORS.log"

        self._logger = logging.getLogger(self.logger_name)
        self._logger.handlers.clear()
        self._logger.propagate = False
        
        self._logger.setLevel(logging.DEBUG if (self.debug_mode or self.test_mode) else logging.INFO)

        if not os.path.exists("logs"):
            os.makedirs("logs")

        ch = logging.StreamHandler()
        fh = logging.FileHandler(logfilename, encoding='utf-8')
        fh_errors = logging.FileHandler(logfilename_errors, encoding='utf-8')

        ch.setLevel(logging.ERROR if self.test_mode else 
                   logging.DEBUG if self.debug_mode else logging.INFO)
        fh.setLevel(logging.DEBUG if (self.debug_mode or self.test_mode) else logging.INFO)
        fh_errors.setLevel(logging.WARNING)

        logging.addLevelName(logging.ERROR, 'Error')
        logging.addLevelName(logging.WARNING, 'Warning')
        logging.addLevelName(logging.INFO, 'Info')
        logging.addLevelName(logging.DEBUG, 'Debug')
        formatter_log = logging.Formatter(
            '%(levelname)-8s | %(filename)-22s | %(funcName)-26s (line %(lineno)-4s) : %(message)s'
        )
        # use custom formatter to get colored output for the command line
        formatter_cmd = CustomFormatter() 
        
        ch.setFormatter(formatter_cmd)
        fh.setFormatter(formatter_log)
        fh_errors.setFormatter(formatter_log)
        for handler in [ch, fh, fh_errors]:
            self._logger.addHandler(handler)

        self._logger.debug("Logger initialized successfully")

    @property
    def logger(self):
        if self._logger is None:
            raise RuntimeError("Logger has not been properly initialized")
        return self._logger

    @classmethod
    def get_logger(cls, name=None):
        if cls._instance is None:
            if name is None:
                raise RuntimeError("Logger not initialized. Create a Logger instance first.")
            return cls(logger_name=name).logger
        return cls._instance
    
    @classmethod
    def change_log_location(cls, new_log_path):
        if cls._instance is None:
            raise RuntimeError("Logger not initialized. Create a Logger instance first.")
        
        instance = cls.get_logger()
        
        os.makedirs(os.path.dirname(new_log_path), exist_ok=True)
        
        file_handlers = [h for h in instance.handlers if isinstance(h, logging.FileHandler)]
        
        for handler in file_handlers:
            log_level = handler.level
            log_formatter = handler.formatter

            instance.removeHandler(handler)
            handler.close()
            
            if not os.path.exists(new_log_path):
                os.makedirs(new_log_path)

            old_log_path = handler.baseFilename
            if os.path.exists(old_log_path):
                shutil.move(old_log_path, new_log_path)
            new_log_file_path = os.path.join(new_log_path, os.path.basename(old_log_path))
            
            new_file_handler = logging.FileHandler(new_log_file_path, encoding='utf-8')
            new_file_handler.setLevel(log_level)
            new_file_handler.setFormatter(log_formatter)

            instance.addHandler(new_file_handler)

        instance.info(f"Log files moved to: {new_log_path}")

    def close_logger(self, logger):
        if logger.hasHandlers():
            logger.handlers.clear()

