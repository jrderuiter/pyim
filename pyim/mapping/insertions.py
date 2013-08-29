
from __future__ import print_function


def write_to_gff(ins_frame):
    with open('ins.gff', 'w') as gff_file:
        gff_str = '{chrom}\t.\tinsertion\t{start}\t{end}\t.\t{strand}\t.\t{info}'
        info_str = 'ID=INS{name}({barcode});NAME={barcode};LP={lp};uLP={ulp}'

        for name, row in ins_frame.iterrows():
            info_str_fmt = info_str.format(name=name, barcode=row['barcode'], lp=row['lp'], ulp=row['unique_lp'])
            gff_str_fmt = gff_str.format(chrom=row['chromosome'],  strand=row['strand'], info=info_str_fmt,
                                     start=row['location'] - 10, end=row['location'] + 10)
            print(gff_str_fmt, file=gff_file)
