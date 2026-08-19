import numpy as np
import pandas as pd
import pytest

from trcss.core import (
    compute_trcss,
    compute_trcss_dataframe,
    fold_interaction,
    validate_inputs,
)


class TestValidateInputs:
    def test_scalar_broadcast(self):
        fd, tx = validate_inputs(0.5, 1.0)
        assert fd.shape == (1,)

    def test_arrays_equal_length(self):
        fd, tx = validate_inputs([-1.0, 0.0, 1.0], [1.0, -1.0, 1.0])
        assert fd.shape == (3,)

    def test_mismatched_length_raises(self):
        with pytest.raises(ValueError):
            validate_inputs([1.0, 2.0], [1.0])

    def test_fd_out_of_range_raises(self):
        with pytest.raises(ValueError):
            validate_inputs([1.5], [1.0])
        with pytest.raises(ValueError):
            validate_inputs([-1.1], [1.0])

    def test_txdir_invalid_raises(self):
        with pytest.raises(ValueError):
            validate_inputs([0.0], [0.0])
        with pytest.raises(ValueError):
            validate_inputs([0.0], [2.0])

    def test_nan_allowed(self):
        fd, tx = validate_inputs([np.nan, 0.5], [1.0, np.nan])
        assert np.isnan(fd[0])
        assert np.isnan(tx[1])


class TestFoldInteraction:
    def test_head_on_maps_to_one(self):
        assert float(fold_interaction(-1.0)) == 1.0

    def test_codirectional_maps_to_zero(self):
        assert float(fold_interaction(+1.0)) == 0.0

    def test_midpoint(self):
        assert float(fold_interaction(0.0)) == pytest.approx(0.5)

    def zero_info_loss_spearman(self):
        signed = np.linspace(-1, 1, 101)
        folded = fold_interaction(signed)
        rho = pd.Series(signed).corr(pd.Series(folded), method="spearman")
        assert rho == pytest.approx(-1.0)


class TestComputeTrcss:
    def test_head_on_plus_strand(self):
        # FD=-1 (leftward), TxDir=+1 (plus gene) -> FD*TxDir=-1 -> TRCSS=1
        assert float(compute_trcss(-1.0, +1.0)) == 1.0

    def test_codirectional_plus_strand(self):
        # FD=+1 (rightward), TxDir=+1 -> +1 -> TRCSS=0
        assert float(compute_trcss(+1.0, +1.0)) == 0.0

    def test_head_on_minus_strand(self):
        # FD=+1 (rightward fork), TxDir=-1 (minus gene, runs left)
        # Gene on minus strand advances to the left. A rightward-moving fork
        # advances to the right.  They advance toward each other -> head-on.
        # FD * TxDir = (+1) * (-1) = -1 -> TRCSS = 1
        assert float(compute_trcss(+1.0, -1.0)) == 1.0

    def test_codirectional_minus_strand(self):
        # FD=-1 (leftward fork), TxDir=-1 (minus gene) -> both move leftward
        # FD * TxDir = +1 -> TRCSS = 0
        assert float(compute_trcss(-1.0, -1.0)) == 0.0

    def test_array_shape(self):
        trcss = compute_trcss([-1.0, 0.0, 1.0], [1.0, 1.0, 1.0])
        assert trcss.shape == (3,)
        np.testing.assert_allclose(trcss, [1.0, 0.5, 0.0])


class TestComputeTrcssDataframe:
    def test_adds_column(self):
        df = pd.DataFrame({
            "fd": [-1.0, 0.0, 1.0],
            "txdir": [1.0, 1.0, 1.0],
        })
        out = compute_trcss_dataframe(df)
        assert "trcss" in out.columns
        np.testing.assert_allclose(out["trcss"].values, [1.0, 0.5, 0.0])

    def test_inplace(self):
        df = pd.DataFrame({"fd": [0.0], "txdir": [-1.0]})
        res = compute_trcss_dataframe(df, inplace=True)
        assert res is None
        assert "trcss" in df.columns
        assert float(df["trcss"].iloc[0]) == pytest.approx(0.5)

    def test_custom_colnames(self):
        df = pd.DataFrame({"fork": [-1.0], "strand": [1.0]})
        out = compute_trcss_dataframe(
            df, fd_col="fork", txdir_col="strand", out_col="TRCSS_final"
        )
        assert "TRCSS_final" in out.columns
        assert float(out["TRCSS_final"].iloc[0]) == 1.0
