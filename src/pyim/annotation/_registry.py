
_registry = {}


def register_annotator(name, aligner):
    _registry[name] = aligner


def get_annotators():
    return dict(_registry)
