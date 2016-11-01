import abc
from pathlib import Path

from pyim.util.path import build_path

from ..external.cutadapt import cutadapt
from .base import Pipeline


class PairedPipeline(Pipeline):
    @abc.abstractclassmethod
    def configure_args(cls, parser):
        parser.add_argument('--reads', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)

    @abc.abstractclassmethod
    def from_args(cls, args):
        raise NotImplementedError()

    @abc.abstractclassmethod
    def run(self, reads_path, work_dir):
        raise NotImplementedError()

    def extract_genomic(self, reads_path, output_base):
        # Ensure output dir exists.
        output_base.parent.mkdir(exist_ok=True)
