__author__ = 'Julian'


class Annotator(object):

    def __init__(self, **kwargs):
        super(Annotator, self).__init__()

    def annotate(self, frame):
        raise NotImplementedError

    @classmethod
    def register_parser(cls, subparsers, name):
        parser = subparsers.add_parser(name=name)

        # Set target class.
        parser.set_defaults(annotator=cls)

        # Add default positional arguments.
        parser.add_argument('input')
        parser.add_argument('output')

        return parser
