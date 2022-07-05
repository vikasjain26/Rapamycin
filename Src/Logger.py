
import logging 

logginFormat = "%(asctime)s- %(levelname)s - %(message)s"

def get_Logger(name):
    return logging.getLogger(name)

def setup_logger(loggetinstance,name,log_file,log_level = logging.INFO):
    """
    """
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(logginFormat)
    
    loggetinstance.addHandler(log_handler)
    loggetinstance.setLevel(log_level)

    #return logger_obj


def close_filehandlers(loginstance):
    """
    """
    h = loginstance.handlers[:]
    for handles in h:
        handles.close()
        loginstance.removeHandler(handles)
    