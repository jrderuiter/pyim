import argparse

from pyim.annotate import get_annotators
from pyim.model import Insertion


def main():
    """Main function for pyim-annotate."""
    args = parse_args()

    insertions = Insertion.from_csv(args.insertions, sep='\t')

    annotator = args.caller.from_args(args)
    annotated = list(annotator.annotate(insertions))

    annotated_frame = Insertion.to_frame(annotated)
    annotated_frame = annotated_frame.sort_values(by='id')
    annotated_frame.to_csv(str(args.output), sep='\t', index=False)


def parse_args():
    """Parses arguments for pyim-annotate."""

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-annotate')
    subparsers = parser.add_subparsers(dest='annotator')
    subparsers.required = True

    # Register pipelines.
    for name, class_ in get_annotators().items():
        cis_parser = subparsers.add_parser(name)
        class_.configure_args(cis_parser)
        cis_parser.set_defaults(caller=class_)

    return parser.parse_args()


if __name__ == '__main__':
    main()
