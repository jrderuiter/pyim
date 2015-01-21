from copy import deepcopy

from pyim.common.alignment.vector import aligners


class SequenceAlignerFactory(object):

    PARALLEL_ALIGNERS = {
        'query': aligners.ParallelQuerySequenceAligner,
        'vector': aligners.ParallelVectorSequenceAligner
    }

    @classmethod
    def make_from_options(cls, opts):
        opts = deepcopy(opts)

        if opts['options'] is not None:
            # Extract options dict from opt for convenience.
            opts_options = opts['options']

            # Replace 'aligner' if present with actual aligner.
            if 'aligner' in opts_options:
                opts_options['aligner'] = cls.make_from_options(opts_options['aligner'])

            # Replace 'aligners' if present with actual aligners.
            if 'aligners' in opts_options:
                opts_options['aligners'] = \
                    [cls.make_from_options(o) for o in opts_options['aligners']]

            # Replace original options dict with instantiated instance.
            opts['options'] = opts_options

        return cls._build_aligner(opts)

    @classmethod
    def _build_aligner(cls, opts):
        class_ = getattr(aligners, opts['type'])
        options = opts['options'] if opts['options'] is not None else {}
        return class_(**options)

    @classmethod
    def make_parallel_from_options(cls, opts, threads, p_type='query'):
        aligner = cls.make_from_options(opts)

        if threads > 1:
            try:
                p_class = cls.PARALLEL_ALIGNERS[p_type]
                aligner = p_class(aligner=aligner, threads=threads)
            except KeyError:
                raise KeyError('Unknown parallel type {}'.format(p_type))

        return aligner
