import argparse

from pyim.alignment.pipelines.sb import SbPipeline


def alignment_main(options):
    if options.pipeline == 'sb':
        SbPipeline.run(options)
    else:
        raise NotImplementedError('Unknown pipeline "%s"' % options.pipeline)


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', required=True)
    ##parser.add_argument('-r', '--reference', required=True)
    parser.add_argument('--pipeline', default='sb')

    return parser.parse_args()

if __name__ == '__main__':
    alignment_main(_parse_args())

