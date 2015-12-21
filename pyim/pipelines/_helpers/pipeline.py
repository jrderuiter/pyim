import collections
import itertools

import skbio
import toolz


@toolz.curry
def print_stats(results, logger=None, header=False):
    print_ = print if logger is None else logger.info

    # Iterate over results, counting statuses.
    status_counts = collections.defaultdict(int)

    for result in results:
        status_counts[result.status.name] += 1
        yield result

    # We're done, so print frequencies!
    if header:
        print_('Extraction stats:')

    total = sum(status_counts.values())
    for status, count in status_counts.items():
        percentage = (count / total) * 100
        print_('{:>18}: {:>8} ({:05.2f}%)'
               .format(status, count, percentage))


@toolz.curry
def write_genomic_sequences(results, file_path, format='fastq',
                            mode='w', **io_kwargs):
    """ Test docstring """
    with skbio.io.open(file_path, mode, **io_kwargs) as file_:
        for result in results:
            skbio.io.write(result.genomic_sequence, into=file_, format=format)
            yield result


@toolz.curry
def build_barcode_map(results, sample_map=None):
    if sample_map is None:
        return {result.genomic_sequence.metadata['id']:
                result.barcode
                for result in results}
    else:
        return {result.genomic_sequence.metadata['id']:
                sample_map[result.barcode]
                for result in results}


def consume(iterator, n=None):
    "Advance the iterator n-steps ahead. If n is none, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # Feed the entire iterator into a zero-length deque
        collections.deque(iterator, maxlen=0)
    else:
        # Advance to the empty slice starting at position n
        next(itertools.islice(iterator, n, n), None)
