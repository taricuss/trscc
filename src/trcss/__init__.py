"""TRCSS: Transcription-Replication Context Score.

A portable, single-feature engineered covariate that captures the joint
geometry of replication fork directionality, transcription strand
orientation, and distance to the nearest transcription start site for
prime editing efficiency prediction.
"""

from .core import (
    compute_trcss,
    compute_trcss_dataframe,
    fold_interaction,
    validate_inputs,
)

__version__ = "1.0.0"
__all__ = [
    "compute_trcss",
    "compute_trcss_dataframe",
    "fold_interaction",
    "validate_inputs",
    "__version__",
]
