"""Core TRCSS computation functions.

TRCSS (Transcription-Replication Context Score) quantifies the local
geometric relationship between replication fork directionality (FD) and
transcription strand orientation (TxDir). The score is defined as:

    TRCSS = (1 - FD * TxDir) / 2

where:
    FD    - Fork Directionality index from OK-seq/Repli-seq, range [-1, +1]
            +1 = fully rightward (Watson-strand Okazaki-dominant)
            -1 = fully leftward (Crick-strand Okazaki-dominant)
    TxDir - Transcription direction, {-1, +1}
            +1 = plus-strand gene
            -1 = minus-strand gene

The folded formulation maps:
    head-on geometry       (FD * TxDir = -1) -> TRCSS = 1
    co-directional geometry (FD * TxDir = +1) -> TRCSS = 0

Empirically, higher TRCSS values are associated with higher prime editing
efficiency in the K562 validation dataset (Mathis et al. 2025).
"""

from __future__ import annotations

from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd


Number = Union[float, int, np.floating, np.integer]


def validate_inputs(
    fd: Union[Number, np.ndarray, pd.Series, Iterable[Number]],
    txdir: Union[Number, np.ndarray, pd.Series, Iterable[Number]],
) -> tuple[np.ndarray, np.ndarray]:
    """Validate and coerce FD and TxDir inputs.

    Parameters
    ----------
    fd : array-like
        Fork directionality values. Should be in [-1, +1].
    txdir : array-like
        Transcription direction values. Should be in {-1, +1}.

    Returns
    -------
    (fd_arr, txdir_arr) : tuple of np.ndarray
        Validated 1-D float arrays of equal length.

    Raises
    ------
    ValueError
        If FD values fall outside [-1, 1] or TxDir values are not in
        {-1, +1}, or if the two arrays have mismatched lengths.
    """
    fd_arr = np.asarray(fd, dtype=np.float64).ravel()
    txdir_arr = np.asarray(txdir, dtype=np.float64).ravel()

    if fd_arr.size != txdir_arr.size:
        raise ValueError(
            f"FD and TxDir arrays must have the same length, "
            f"got {fd_arr.size} vs {txdir_arr.size}."
        )

    fd_mask = ~(np.isnan(fd_arr) | (fd_arr >= -1.0) & (fd_arr <= 1.0))
    if fd_arr.size and np.any(fd_mask & ~np.isnan(fd_arr)):
        bad = fd_arr[~np.isnan(fd_arr) & fd_mask]
        if bad.size:
            raise ValueError(
                f"FD values must be in [-1, +1] or NaN. "
                f"Found out-of-range: {bad[:5]}."
            )

    txdir_valid = np.isnan(txdir_arr) | np.isclose(txdir_arr, 1.0) | np.isclose(txdir_arr, -1.0)
    if txdir_arr.size and not np.all(txdir_valid):
        bad = txdir_arr[~txdir_valid]
        raise ValueError(
            f"TxDir values must be -1, +1, or NaN. "
            f"Found invalid: {bad[:5]}."
        )

    return fd_arr, txdir_arr


def fold_interaction(
    fd_times_txdir: Union[Number, np.ndarray, pd.Series, Iterable[Number]],
) -> np.ndarray:
    """Rescale the signed FD * TxDir product to the [0, 1] interval.

    This is the geometrically correct folding step with zero information
    loss: the Spearman correlation between the signed interaction and the
    folded TRCSS is exactly -1.0.

    Parameters
    ----------
    fd_times_txdir : array-like
        Pointwise product of FD and TxDir. Range is [-1, +1].

    Returns
    -------
    trcss : np.ndarray
        Folded TRCSS values in [0, 1], same shape as input.
    """
    arr = np.asarray(fd_times_txdir, dtype=np.float64)
    return (1.0 - arr) / 2.0


def compute_trcss(
    fd: Union[Number, np.ndarray, pd.Series, Iterable[Number]],
    txdir: Union[Number, np.ndarray, pd.Series, Iterable[Number]],
) -> np.ndarray:
    """Compute TRCSS from FD and TxDir.

    Parameters
    ----------
    fd : array-like
        Fork directionality values in [-1, +1].
    txdir : array-like
        Transcription direction values in {-1, +1}.

    Returns
    -------
    trcss : np.ndarray
        TRCSS values in [0, 1]. Same length as inputs.

    Examples
    --------
    Perfect head-on geometry:
    >>> compute_trcss(-1.0, +1.0)
    array([1.])

    Perfect co-directional geometry:
    >>> compute_trcss(+1.0, +1.0)
    array([0.])
    """
    fd_arr, txdir_arr = validate_inputs(fd, txdir)
    return fold_interaction(fd_arr * txdir_arr)


def compute_trcss_dataframe(
    df: pd.DataFrame,
    fd_col: str = "fd",
    txdir_col: str = "txdir",
    out_col: str = "trcss",
    inplace: bool = False,
) -> Optional[pd.DataFrame]:
    """Append a TRCSS column to a pandas DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Input data.
    fd_col : str, default "fd"
        Name of the column containing fork directionality values.
    txdir_col : str, default "txdir"
        Name of the column containing transcription direction values.
    out_col : str, default "trcss"
        Name of the output column.
    inplace : bool, default False
        If True, modify `df` in place. Otherwise return a copy.

    Returns
    -------
    pd.DataFrame or None
        DataFrame with the new TRCSS column, or None if ``inplace=True``.
    """
    if not inplace:
        df = df.copy()
    trcss = compute_trcss(df[fd_col].values, df[txdir_col].values)
    df[out_col] = trcss
    if inplace:
        return None
    return df
