import abc


class CisCaller(abc.ABC):
    """Base CisCaller class."""

    def __init__(self):
        pass

    @abc.abstractmethod
    def call(self, insertions):
        """Calls CIS sites for insertions.

        Returns a tuple of a CisSet and a Dataframe with two columns
        describing the mapping of insertions -> cis sites.
        """

        raise NotImplementedError()
