from collections import namedtuple

import pandas


Sequence = namedtuple('Sequence', ['name', 'sequence'])

class Insertion(object):
    """docstring for Insertion"""
    
    def __init__(self, name, seqname, location, strand, barcode, 
                 sample=None, metadata=None):
        super(Insertion, self).__init__()

        if metadata is None:
            metadata = {}
        
        self.name = name
        self.seqname = seqname
        self.location = location
        self.strand = strand
        self.barcode = barcode
        self.sample = sample
        self.metadata = metadata
    
    def __repr__(self):
        base_str = ("Insertion(name='{name}', seqname='{seqname}', "
                    "location={location}, strand={strand}, "
                    "barcode='{barcode}', sample='{sample}', "
                    "metadata={metadata})")
        return base_str.format(**self.__dict__) 

    @property
    def _dict(self):
        dict_ = self.metadata
        dict_.update(self.__dict__)
        del dict_['metadata']
        return dict_

    @classmethod
    def to_frame(Class, insertions):
        ins_frame = pandas.DataFrame((ins._dict for ins in insertions))

        ## Order our columns accordingly.
        ordered_cols = ['name', 'seqname', 'location', 'strand', 'barcode', 'sample']
        ordered_cols += [c for c in ins_frame.columns if c not in ordered_cols]

        return ins_frame[ordered_cols] 

    @classmethod
    def to_file(Class, insertions, file_path, sep='\t', sort_cols=['seqname', 'location']):
        ## Convert insertions to a data frame.
        ins_frame = Class.to_frame(insertions)

        ## Sort if required.
        if sort_cols is not None:
            ins_frame.sort(columns=sort_cols, inplace=True)
    
        ## Write final output to frame.
        ins_frame.to_csv(file_path, sep=sep, header=True, index=False)
