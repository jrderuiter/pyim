
import pandas


def map_insertions_to_samples(insertions, barcode_mapping_file, drop_barcode=True):
    barcode_map_frame = pandas.read_csv(barcode_mapping_file)
    barcode_map = dict(zip(barcode_map_frame['Barcode'], barcode_map_frame['Sample']))
    insertions['sample'] = insertions['barcode'].map(barcode_map)

    if drop_barcode:
        del insertions['barcode']

    return insertions