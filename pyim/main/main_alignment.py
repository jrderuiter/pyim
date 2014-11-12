import argparse
import yaml
from pyim.alignment.pipelines.sb import SbPipeline

PIPELINES = {
    'sb': SbPipeline
}


def alignment_main(options):
    ## Look-up the requested pipeline.
    try:
        pipeline = PIPELINES[options.pipeline]
    except KeyError:
        raise KeyError('Unknown pipeline "{}"'.format(options.pipeline))

    ## Load pipeline config.
    with open(options.config, 'r') as config_file:
        config = yaml.load(config_file)

    ## Actually run the pipeline!
    pipeline.run(options.input, options.output, config)


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)

    parser.add_argument('-p', '--pipeline', default='sb')
    parser.add_argument('-c', '--config', required=True)

    return parser.parse_args()


if __name__ == '__main__':
    alignment_main(_parse_args())

