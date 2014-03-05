
import argparse
import pandas


def main():
    args = _parse_args()

    if not len(args.insertion_sets) == len(args.run_names):
        raise ValueError("The number of names is not equal to the number of sets!")

    # Load and merge frames, injecting run names as extra column
    ins_frames = []
    for run_name, set_file in zip(args.set_names, args.insertion_sets):
        ins_frame = pandas.read_csv(set_file, sep='\t')
        ins_frame['run'] = run_name
        ins_frames.append(ins_frame)
    ins_merged = pandas.concat(ins_frames, ignore_index=True)

    # Mask samples if requested. Argument 'keep_only_masked'
    # determines if samples in the mask are kept or discarded.
    ins_in_mask = ins_merged['sample'].isin(args.mask)
    if args.keep_only_masked:
        ins_merged = ins_merged[ins_in_mask]
    else:
        ins_merged = ins_merged[~ins_in_mask]

    # Reset the ids in the merged frame to ensure ids are unique
    ins_merged['id'] = ['INS_%d' % d for d in range(1, len(ins_merged)+1)]

    ins_merged.to_csv(args.output, sep='\t', index=False)


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', dest='insertion_sets', required=True, nargs='+')
    parser.add_argument('-n', '--names', dest='run_names', required=True, nargs='+')
    parser.add_argument('-o', '--output', dest='merged_output', required=True)
    parser.add_argument('-m', '--mask', dest='masked_files', default=[], nargs='+')
    parser.add_argument('--keep-mask', dest='keep_only_masked', default=False, action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    main()



#def setup_parser(root=None):
    #if root is not None:
    #    parser = root.add_parser('merge')
    #else:
    #    parser = argparse.ArgumentParser()
