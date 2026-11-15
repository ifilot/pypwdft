from importlib.metadata import PackageNotFoundError as _PackageNotFoundError
from importlib.metadata import version as _distribution_version

from .api import (
    BasisInfo,
    DFTResult,
    EnergyComponents,
    GTH,
    OrbitalSet,
    PWDFT,
    SCFInfo,
    SCFSettings,
    Structure,
)

__all__ = [
    "BasisInfo",
    "DFTResult",
    "EnergyComponents",
    "GTH",
    "OrbitalSet",
    "PWDFT",
    "SCFInfo",
    "SCFSettings",
    "Structure",
    "__version__",
]

try:
    __version__ = _distribution_version("pypwdft")
except _PackageNotFoundError:
    # The distribution metadata is unavailable when importing an uninstalled
    # source tree. Installed and editable packages obtain this from
    # ``pyproject.toml`` through their generated metadata.
    __version__ = "0+unknown"
