
def aligner_from_options(options):
    from . import aligners
    class_ = getattr(aligners, options['type'])
    return class_(**options['options'])
