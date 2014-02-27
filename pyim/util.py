
import logging


def setup_logging():
    root = logging.getLogger()
    if root.handlers:
        for handler in root.handlers:
            root.removeHandler(handler)
    logging.basicConfig(format='%(asctime)s %(name)s \t %(message)s', level=logging.INFO)


def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

