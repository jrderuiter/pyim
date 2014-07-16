import argparse
import pandas


def main(options):

    # Load and merge frames, injecting run names as extra column
    ins_frames = []
    for ins_file in options.insertion_sets:
        ins_frame = pandas.read_csv(ins_file, sep='\t')
        ins_frames.append(ins_frame)
    ins_merged = pandas.concat(ins_frames, ignore_index=True)

    ## First we check if all masked samples are indeed present
    ## in the data set to avoid issues due to spelling errors etc.
    sample_set = set(ins_merged['sample'])
    for sample in options.masked_samples:
        if sample not in sample_set:
            print("WARNING: Requested masked sample %s was not encountered in the dataset." % sample)

    ## Mask samples if requested. Argument 'keep_only_masked'
    ## determines if samples in the mask are kept or discarded.
    ins_in_mask = ins_merged['sample'].isin(options.masked_samples)
    if options.keep_only_masked:
        ins_merged = ins_merged[ins_in_mask]
    else:
        ins_merged = ins_merged[~ins_in_mask]

    # Reset the ids in the merged frame to ensure ids are unique
    if any(ins_merged['id'].duplicated()):
        raise ValueError('Resulting merged frame contains duplicate ids.')

    ins_merged.to_csv(options.merged_output, sep='\t', index=False)


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--insertion-sets', dest='insertion_sets', required=True, nargs='+')
    parser.add_argument('-o', '--output-set', dest='merged_output', required=True)
    parser.add_argument('-m', '--masked-samples', dest='masked_samples', default=[], nargs='+')
    parser.add_argument('--masked-only', dest='keep_only_masked', default=False, action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    main(_parse_args())



#def setup_parser(root=None):
    #if root is not None:
    #    parser = root.add_parser('merge')
    #else:
    #    parser = argparse.ArgumentParser()
