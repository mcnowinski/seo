import logging
import logging.handlers
import time

def get_logger(name):   
    #configure logging
    log_file='./' + name
    logger = logging.getLogger(name)
    logging.Formatter.converter = time.gmtime
    logger.setLevel(logging.DEBUG)
    handler = logging.handlers.RotatingFileHandler(log_file, maxBytes=5*1024*1024, backupCount=10)
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s\t%(message)s'))
    logger.addHandler(handler)

    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s\t%(message)s'))
    logger.addHandler(handler)

    return logger

def get_dumper(name):   
    #configure logging
    dump_file='./' + name
    dumper = logging.getLogger(name)
    logging.Formatter.converter = time.gmtime
    dumper.setLevel(logging.DEBUG)
    handler = logging.handlers.RotatingFileHandler(dump_file, maxBytes=5*1024*1024, backupCount=10)
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s\t%(message)s'))
    dumper.addHandler(handler)

    return dumper