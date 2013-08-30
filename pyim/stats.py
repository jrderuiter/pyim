__author__ = 'j.d.ruiter'

import os
import pandas

from matplotlib import pyplot as plt


def hist_alignment_type(alignments, unmapped=None, name='', filepath=None):
    aln_counts = alignments['type'].value_counts()

    if unmapped is not None:
        aln_counts = aln_counts.append(pandas.Series(len(unmapped), index=['unmapped']))

    aln_counts.plot(kind='bar')

    title_str = "%s Alignment types" % name if name else "Alignment types"
    plt.title(title_str); plt.xlabel('Alignment type'); plt.ylabel('Count')
    _show(filepath)

    return aln_counts


def hist_barcode_alignment(bc_alignments, fasta_seqs=None, filepath=None):
    bc_alignments = pandas.concat(bc_alignments.values(), ignore_index=True)

    bc_counts = bc_alignments['target_name'].value_counts()
    bc_counts = bc_counts.iloc[bc_counts.index.argsort()]

    if fasta_seqs is not None:
        read_is_bc_mapped = pandas.Series([s.name for s in fasta_seqs]).isin(bc_alignments['query_name'])
        bc_unmapped_num = len((~read_is_bc_mapped).nonzero()[0])
        bc_counts = bc_counts.append(pandas.Series(bc_unmapped_num, index=['unmapped']))

    bc_counts.plot(kind='bar')
    _show(filepath)

    return bc_counts

def _show(filepath):
    if filepath is None:
        plt.show()
    else:
        save(filepath)


def save(path, close=True, verbose=False):
    """Save a figure from pyplot.

    Parameters
    ----------
    path : string
        The path (and filename, without the extension) to save the
        figure to.

    ext : string (default='pdf')
        The file extension. This must be supported by the active
        matplotlib backend (see matplotlib.backends module).  Most
        backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.

    close : boolean (default=True)
        Whether to close the figure after saving.  If you want to save
        the figure multiple times (e.g., to multiple formats), you
        should NOT close it in between saves or you will have to
        re-plot it.

    verbose : boolean (default=False)
        Whether to print information about when and where the image
        has been saved.

    """

    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = os.path.split(path)[1]

    if directory == '':
        directory = '.'

    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)

    # The final path to save to
    save_path = os.path.join(directory, filename)

    if verbose:
        print("Saving figure to '%s'..." % save_path),

    # Actually save the figure
    plt.savefig(save_path)

    # Close it
    if close:
        plt.close()