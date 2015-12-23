import logging
import pkg_resources


logging.basicConfig(
        format='%(asctime)-15s %(message)s',
        datefmt='[%Y-%m-%d %H:%M:%S]',
        level=logging.INFO)


def print_header(logger):
    version = pkg_resources.require('pyim')[0].version
    header_str = ' PyIM ({}) '.format(version)
    logger.info('{:-^40}'.format(header_str))


def print_footer(logger):
    logger.info('{:-^40}'.format(' Done! '))
