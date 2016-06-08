import logging
import pkg_resources


logging.basicConfig(
        format='%(asctime)-15s %(levelname)-10s %(message)s',
        datefmt='[%Y-%m-%d %H:%M:%S]',
        level=logging.INFO)


def print_header(logger, command=None):
    version = pkg_resources.require('pyim')[0].version

    if command is None:
        header_str = ' PyIM ({}) '.format(version)
    else:
        header_str = ' PyIM {} ({}) '.format(command, version)

    logger.info('{:-^60}'.format(header_str))


def print_footer(logger):
    logger.info('{:-^60}'.format(' Done! '))
