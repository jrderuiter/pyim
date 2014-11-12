from itertools import repeat, chain


def chunks(l, n):
    """Divide a list l into chunks of size n.

    :Parameters:
      - `l` (list) - List of objects to divide into chunks.
      - `n` (int)  - Chunk size.

    :Returns:
        List l divided into chunks of at most size n.

    :Returns Type:
        List of lists.

    :Examples:
        chunked_list = chunks([1, 2, 3], 2)

    """

    return [l[i:i+n] for i in range(0, len(l), n)]


def repeat_list(l, n):
    return chain(*[repeat(el, num) for el, num in zip(l, n)])


def flatten_list(l):
    return [item for sub_list in l for item in sub_list]