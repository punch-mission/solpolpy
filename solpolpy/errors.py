class SolpolpyError(Exception):
    pass


class TooFewFilesError(SolpolpyError):
    pass


class UnsupportedInstrumentError(SolpolpyError):
    pass


class UnsupportedTransformationError(SolpolpyError):
    pass


class MissingAlphaError(SolpolpyError):
    pass


class InvalidDataError(SolpolpyError):
    pass
