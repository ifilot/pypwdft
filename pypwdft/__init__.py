from importlib.metadata import PackageNotFoundError, version

from .pypwdft import PyPWDFT
from .psystem import PeriodicSystem
from .system_builder import SystemBuilder
from .gth import GTHPseudopotential

try:
    __version__ = version("pypwdft")
except PackageNotFoundError:
    # The distribution metadata is unavailable when importing an uninstalled
    # source tree. Installed and editable packages obtain this from
    # ``pyproject.toml`` through their generated metadata.
    __version__ = "0+unknown"
