import abc
from pathlib import Path

_registry = {}


def register_pipeline(name, pipeline):
    _registry[name] = pipeline


def get_pipelines():
    return dict(_registry)


class Pipeline(abc.ABC):
    def __init__(self):
        pass

    @abc.abstractclassmethod
    def configure_args(cls, parser):
        parser.add_argument('--reads', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)

    @abc.abstractclassmethod
    def extract_args(cls, args):
        raise NotImplementedError()

    @classmethod
    def from_args(cls, args):
        return cls(**cls.extract_args(args))

    @abc.abstractclassmethod
    def run(self, reads_path, work_dir):
        raise NotImplementedError()
