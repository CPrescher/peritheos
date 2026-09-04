"""Shock-path EOS API under the common :mod:`peritheos.eos` namespace.

The top-level :mod:`peritheos.hugoniot` module remains the compatibility import
path for the initial release of this API.
"""

from peritheos.hugoniot import HugoniotBase, HugoniotState, LinearUsUpHugoniot

__all__ = ["HugoniotBase", "HugoniotState", "LinearUsUpHugoniot"]
