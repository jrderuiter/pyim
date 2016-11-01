

from ..common.cutadapt import cutadapt

def shear_splink(reads, barcodes):

    # De-multiplex
    sample_files = _demultiplex(reads, barcodes)

    # Filter for contaminants


    # Select for SB and T7



def _demultiplex(reads, output, barcodes):
    options = {
        '-g': ('file:' + str(barcodes), )
        '--discard-untrimmed': ()
    }

    cutadapt(reads_path, output_path, options=options)


def _extract_genomic(reads, transposon, linker, contaminants):

    # Filter for contaminants.
    options = {
        '-g': ('file:' + str(contaminants), )
        '--discard-trimmed': ()
    }

    cutadapt(reads_path, tmp_path, options)

    # Select for and remove transposon sequence.

    # Select for and remove linker sequence.
