import logging

logger = logging.getLogger('sldb')

def setup_default_logging(level=logging.INFO):
    logger.setLevel(level)
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'
        )
    )
    logger.addHandler(handler)
