try:
    from ._version import version
    __version__ = version
except ImportError:
    from importlib.metadata import version, PackageNotFoundError

    try:
        __version__ = version("woldrnaseq")
    except PackageNotFoundError:
        # package is not installed
        pass
