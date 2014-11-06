from os import path

import rpy2.robjects as robjects
import pandas.rpy.common as com
import pandas

KCRBM_DIR = path.join(path.dirname(__file__), '../../../lib/kcrbm')


def kcrbm_assign_genes(insertions):
    kcrbm_ins_frame = insertions_to_kcrbm_frame(insertions)
    kcrbm_assignment = run_kcrbm(kcrbm_ins_frame)
    return reformat_kcrbm_result(kcrbm_assignment)


def insertions_to_kcrbm_frame(insertions):
    # KCRBM idata format:
    #		Dataframe:
    #	 		$chr (numeric): chromosome number, with X and Y indicated by numbers
    #		  	$base (numeric): insertion locus (bp)
    #		  	$ori (numeric): insertion orientation, sense or antisense (+1/-1)
    #		NB if idata has a column "id", then it has to be a unique identifier for
    #		an insertion. If no such column is present, a column "id" will be added.

    kcrbm_frame = insertions[['id', 'chromosome', 'location', 'strand']].copy()
    kcrbm_frame.columns = ['id', 'chr', 'base', 'ori']

    # Map columns to values expected by KCRBM. Chromosomes are expected to be
    # encoded as integers, with X and Y corresponding to integer values of 20
    # and 21 respectively. The orientation of the transposon is encoded as -1
    # or +1, corresponding to the forward and reverse strand respectively.
    chr_map = dict(zip(map(str, range(1, 19+1)) + ['X', 'Y'], range(1, 21+1)))
    kcrbm_frame['chr'] = kcrbm_frame['chr'].map(chr_map)
    kcrbm_frame['ori'] = kcrbm_frame['ori'].map({'+': 1, '-': -1})

    if pandas.isnull(kcrbm_frame['chr']).any():
        raise ValueError('Encountered unknown chromosomes in the insertion frame!')

    return kcrbm_frame


def run_kcrbm(insertions_kcrbm):
    r_frame = com.convert_to_r_dataframe(insertions_kcrbm)
    robjects.globalenv['idata'] = r_frame

    robjects.r("""
        source('{base_dir}/r/kcrbm.r')
        source('{base_dir}/r/kcrbm_tools.r')
        load('{base_dir}/data/edata.mm10.Rdata')

        idata$ins_id <- idata$id    # Ensure id is copied to output

        t1 <- Sys.time()
        iset.gene <- kcrbm(
            edata = edata,
            idata = idata,
            rules = "SB",           # Parameter values used for SB in NAR paper
            map_to = "genes"
        )
        t2 <- Sys.time()
        print(difftime(t2,t1))

        # Lookup single gene CTGs and record in matrix
        iset.ctg <- selectSingleTargets(iset.gene, type='ctg')

        iset.gene$single_ctg <- FALSE
        iset.gene[row.names(iset.ctg), 'single_ctg'] <- TRUE

        #iset.nearest <- selectSingleTargets(iset.gene, type='nearest')
    """.format(base_dir=KCRBM_DIR))

    iset_gene = com.load_data('iset.gene')

    return iset_gene


def reformat_kcrbm_result(kcrbm_assignment):
    sel_columns = ['ins_id', 'clusterpeak', 'clusterorientation', 'ids', 'n_insertions',
                   'ensid', 'd2gss', 'mechanism', 'd2gts', 'single_ctg']
    renamed_columns = ['id', 'cluster_peak', 'cluster_ori', 'cluster_ids', 'cluster_size',
                       'gene_id', 'd2gss', 'mechanism', 'd2gts', 'single_ctg']

    frame = kcrbm_assignment[sel_columns].copy()
    frame.columns = renamed_columns

    return frame



