import argparse
import logging

from pyim.align.pipelines import get_pipelines

logging.basicConfig(
    format='[%(asctime)-15s]  %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def main():
    """Main function for pyim-align."""

    args = parse_args()

    # Run pipeline.
    reads2 = args.reads2 if hasattr(args, 'reads2') else None

    pipeline = args.pipeline.from_args(args)
    pipeline.run(read_path=args.reads,
                 output_dir=args.output_dir,
                 read2_path=reads2)


def parse_args():
    """Parses arguments for pyim-align."""

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-align')
    subparsers = parser.add_subparsers(dest='pipeline')
    subparsers.required = True

    # Register pipelines.
    pipelines = get_pipelines()

    for pipeline_name in sorted(pipelines.keys()):
        pipeline_class = pipelines[pipeline_name]

        pipeline_parser = subparsers.add_parser(pipeline_name)
        pipeline_class.configure_args(pipeline_parser)
        pipeline_parser.set_defaults(pipeline=pipeline_class)

    return parser.parse_args()


if __name__ == '__main__':
    main()
