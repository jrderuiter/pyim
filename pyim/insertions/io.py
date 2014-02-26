
from __future__ import print_function

def write_to_gff(insertions, file_path):
    id_field = 'sample' if 'sample' in insertions.columns else 'barcode'

    with open(file_path, 'w') as gff_file:
        gff_str = '{chrom}\t.\tinsertion\t{start}\t{end}\t.\t{strand}\t.\t{info}'
        info_str = 'ID={name}({id});NAME={id};LP={lp};uLP={ulp}'

        for i, (_, row) in enumerate(insertions.iterrows()):
            name = row['ins_id'] if 'ins_id' in row else 'INS_' + str(i)

            info_str_fmt = info_str.format(name=name, id=row[id_field], lp=row['lp'], ulp=row['unique_lp'])
            gff_str_fmt = gff_str.format(chrom=row['chromosome'],  strand=row['strand'], info=info_str_fmt,
                                     start=row['location'] - 10, end=row['location'] + 10)

            print(gff_str_fmt, file=gff_file)

    return file_path
