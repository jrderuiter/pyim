import argparse

from frozendict import frozendict
import toolz

from pyim.cis import get_callers
from pyim.model import Insertion, CisSite


def main():
    """Main function for pyim-cis."""

    args = parse_args()
    caller = args.caller.from_args(args)

    # Identify CIS sites.
    insertions = list(Insertion.from_csv(args.insertions, sep='\t'))
    cis_sites, cis_mapping = caller.call(insertions=insertions)

    # Annotate insertions.
    annotated_ins = _annotate_insertions(insertions, cis_mapping)

    # Write outputs.
    Insertion.to_csv(args.output, annotated_ins, sep='\t', index=False)

    if args.output_sites is None:
        cis_path = args.output.with_suffix('.sites.txt')
    else:
        cis_path = args.output_sites

    CisSite.to_csv(cis_path, cis_sites, sep='\t', index=False)


def _annotate_insertions(insertions, cis_map):
    """Annotates insertions with CIS sites using given mapping."""

    for insertion in insertions:
        cis_ids = cis_map.get(insertion.id, set())

        for cis_id in cis_ids:
            cis_metadata = {'cis_id': cis_id}
            new_metadata = toolz.merge(insertion.metadata, cis_metadata)
            yield insertion._replace(metadata=frozendict(new_metadata))


def parse_args():
    """Parses arguments for pyim-cis."""

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-cis')
    subparsers = parser.add_subparsers(dest='caller')
    subparsers.required = True

    # Register pipelines.
    for name, class_ in get_callers().items():
        cis_parser = subparsers.add_parser(name)
        class_.configure_args(cis_parser)
        cis_parser.set_defaults(caller=class_)

    return parser.parse_args()


if __name__ == '__main__':
    main()
