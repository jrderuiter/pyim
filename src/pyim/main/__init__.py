import abc


class Command(abc.ABC):
    """Base command class."""

    name = 'base'

    @abc.abstractmethod
    def configure(self, parser):
        """Configures the argument parser for the command."""
        pass

    @abc.abstractmethod
    def run(self, args):
        """Runs the command."""
        pass

    @classmethod
    def available_commands(cls):
        """Returns dict of available commands."""

        return {
            class_.name: class_()
            for class_ in cls._get_subclasses() if class_.name != 'base'
        }

    @classmethod
    def _get_subclasses(cls):
        for subclass in cls.__subclasses__():
            yield from subclass._get_subclasses()
            yield subclass
