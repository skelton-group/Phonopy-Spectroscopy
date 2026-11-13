# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Class for storing and working with dielectric functions."""

# -------
# Imports
# -------

import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

from ..constants import ZERO_TOLERANCE

from ..units import (
    get_supported_frequency_units,
    get_frequency_unit_text_label,
    get_frequency_unit_plot_label,
)

from ..utility.diagonalisation import reorder_with_branch_tracking

# ----------------------
# ResponseFunction class
# ----------------------


class ResponseFunction:
    """Class for storing and working with general energy-dependent response
    functions."""

    def __init__(self, e, response_func, e_units="ev"):
        """Create a new instance of the `ResponseFunction` class.

        Parameters
        ----------
        e : float or array_like, optional
            Energy in `e_units` (array shape: `(O,)`).
        response_func : array_like
            Response function (shape: `(O, D, D)`, `(D, D)` or `(O,)`).
        e_units : str, optional
            Units of `e` (default: "ev").
        """

        e = np.atleast_1d(np.array(e, dtype=np.float64, copy=True))

        if e.ndim != 1:
            raise ValueError(
                "e must be a scalar or an array_like with shape (O,)."
            )

        response_func = np.atleast_1d(
            np.array(response_func, dtype=np.complex128, copy=True)
        )

        match response_func.ndim:
            case 3:
                # Assume (O, D, D).
                pass
            case 2:
                # Assume (D, D) -> (O, D, D).
                response_func = response_func[np.newaxis, :, :]
            case 1:
                # Assume (O,) -> (O, D, D).
                response_func = response_func[:, np.newaxis, np.newaxis]
            case _:
                raise ValueError(
                    "response_func must be an array_like with shape "
                    "(O, D, D), (D, D) or (O,)."
                )

        if response_func.shape[-2] != response_func.shape[-1]:
            raise ValueError(
                "The outermost two dimensions of eps after reshaping "
                "must match."
            )

        e_units = str(e_units).lower()

        if e_units not in get_supported_frequency_units():
            raise ValueError(f'e_units="{e_units}" are not supported.')

        e.flags.writeable = False
        response_func.flags.writeable = False

        self._e = e
        self._response_func = response_func

        self._e_units = e_units

    @property
    def energy(self):
        """numpy.ndarray : Energies in `energy_units` (shape: `(O,)`)."""
        return self._e

    @property
    def response_function(self):
        """numpy.ndarray : Response function (shape: `(O, D, D)`)."""
        return self._response_func

    @property
    def energy_units(self):
        """str : Units of `energy`."""
        return self._e_units

    @property
    def num_dims(self):
        """int : Number of dimensions D."""
        return self._response_func.shape[-1]

    @property
    def energy_unit_text_label(self):
        """str : Energy unit label suitable for plain-text output."""
        return get_frequency_unit_text_label(self._e_units)

    @property
    def energy_unit_plot_label(self):
        """str : Energy unit label suitable for plotting (may contain
        TeX strings)."""
        return get_frequency_unit_plot_label(self._e_units)

    def eigenmodes(self, branch_tracking=True):
        """Diagonalise the response function and return the eigenmodes
        with optional branch tracking.

        Parameters
        ----------
        branch_tracking : bool, optional
            If `True`, reorder the eigenmodes by performing
            "branch tracking" over energies (default: `True`).

        Returns
        -------
        ems : tuple of numpy.ndarray
            `(evals, evecs)` tuple (shapes: `(O, D)`, `(O, D, D)`).
        """

        if self.num_dims == 1:
            # 1D dielectric function.

            return (
                np.copy(self._response_func[:, 0, 0]),
                np.ones((len(self._e), 1, 1), dtype=np.float64),
            )

        evals, evecs = np.linalg.eig(self._response_func)

        # If branch_tracking is set, reorder evals/evecs.

        if branch_tracking:
            evals, evecs = reorder_with_branch_tracking(evals, evecs)

        evecs = np.real_if_close(evecs, tol=ZERO_TOLERANCE)

        return (evals, evecs)

    @classmethod
    def align_energies(cls, response_funcs, e=None, **kwargs):
        """Align the energies of a set of response functions.

        Parameters
        ----------
        response_funcs : array_like
            Response function(s) to align.
        e : array_like or None, optional
            Energies to align to (default: automatically determined).
        **kwargs : any
            Optional keyword arguments to cls.

        Returns
        -------
        aligned : array_like
            Response functions with aligned energies (same type as the
            calling class).

        Notes
        -----
        Where requiured, linear interpolation is performed using the
        SciPy `interp1d()` routine. If set, all energies in `e` must
        fall within the bounds of all the response functions in
        `response_funcs` (i.e. no extrapolation is performed). If `e` is
        not set, the function will attempt to find the smallest common
        range of energies and smallest step to use for interpolation.

        The returned objects in `aligned` are instances of the calling
        class.
        """

        ret_scalar = np.ndim(response_funcs) == 0

        response_funcs = np.atleast_1d(
            np.asanyarray(response_funcs, dtype=object)
        )

        num_dims = None

        for idx in np.ndindex(*response_funcs.shape):
            if num_dims is None:
                num_dims = response_funcs[idx].num_dims

            if response_funcs[idx].num_dims != num_dims:
                raise ValueError(
                    "All response functions must have the same number "
                    "of dimensions D."
                )

        e_units = None

        for idx in np.ndindex(*response_funcs.shape):
            if e_units is None:
                e_units = response_funcs[idx].energy_units

            if response_funcs[idx].energy_units != e_units:
                # We could handle this, but it almost certainly
                # indicates a mistake.

                raise ValueError(
                    "All response functions must have the same energy_units."
                )

        if e is None:
            # Default to matching the "common" range of energies with
            # the smallest step size.

            e_min = max(f.energy.min() for f in response_funcs.flat)
            e_max = min(f.energy.max() for f in response_funcs.flat)

            if np.isclose(e_min, e_max, atol=ZERO_TOLERANCE):
                # Single energy.

                e = np.array([(e_min + e_max) / 2.0], dtype=np.float64)
            else:
                e_step = min(
                    (f.energy[1:] - f.energy[:-1]).min()
                    for f in response_funcs.flat
                )

                e = np.arange(e_min, e_max + e_step / 10.0, e_step)

            if len(e) == 0:
                raise RuntimeError(
                    "Failed to automatically determine a common set of "
                    "energies for interpolation - one or more response "
                    "functions may not be energy dependent, the "
                    "ranges/steps may not overlap, or e may need to be "
                    "manually specified."
                )

        # Generate a new set of DielectricFunction objects with the
        # energies interpolated to the new axis.

        response_funcs_new = np.zeros_like(response_funcs, dtype=object)

        for idx in np.ndindex(*response_funcs.shape):
            f = response_funcs[idx]

            if len(f.energy) == len(e) and np.allclose(f.energy, e):
                # Energies already aligned.

                if type(f) == cls:
                    response_funcs_new[idx] = f
                else:
                    response_funcs_new[idx] = cls(
                        e,
                        f.response_function,
                        e_units=f.energy_units,
                        **kwargs,
                    )
            else:
                # Interpolation required.

                interp_func = interp1d(
                    f.energy, f.response_function, kind="linear", axis=0
                )

                response_funcs_new[idx] = cls(
                    e, interp_func(e), e_units=f.energy_units, **kwargs
                )

        return response_funcs_new[0] if ret_scalar else response_funcs_new


# ------------------------
# DielectricFunction class
# ------------------------


class DielectricFunction(ResponseFunction):
    """Class for storing and working with energy-dependent 1D, 2D or 3D
    dielectric functions."""

    epsilon_unit_text_label = r"\eps / \eps_0"
    r"""Units of \epsilon suitable for plain-text output."""

    epsilon_unit_plot_label = r"$\varepsilon$ / $\varepsilon_0$"
    r"""Units of \epsilon suitable for plotting (contains TeX strings)."""

    def __init__(self, e, eps, e_units="ev"):
        """Create a new instance of the `DielectricFunction` class.

        Parameters
        ----------
        e : float or array_like, optional
            Energy in `e_units` (array shape: `(O,)`).
        eps : array_like
            Dielectric function in units of relative permittivity
            (shape: `(O, D, D)`, `(D, D)` or `(O,)`).
        e_units : str, optional
            Units of `e` (default: "ev").
        """

        super().__init__(e, eps, e_units=e_units)

        if self.num_dims > 3:
            raise ValueError("D must be 1, 2 or 3.")

    @property
    def epsilon(self):
        """numpy.ndarray : Dielectric function in units of relative
        permittivity (shape: `(O, D, D)`)."""
        return self._response_func

    def dielectric_function(self):
        """Return the dielectric function as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the spectrum.
        """

        d = {"energy": self.energy}

        if self.num_dims == 1:
            # Scalar dielectric function.

            d[f"epsilon_re"] = self.epsilon[:, 0, 0].real
            d[f"epsilon_im"] = self.epsilon[:, 0, 0].imag
        else:
            dim_labels = "xyz"[: self.num_dims]

            for idx, x in enumerate(dim_labels):
                for idy, y in enumerate(dim_labels):
                    eps_re = self.epsilon[:, idx, idy].real
                    d[f"epsilon_re_{x}{y}"] = eps_re

            for idx, x in enumerate(dim_labels):
                for idy, y in enumerate(dim_labels):
                    eps_im = self.epsilon[:, idx, idy].imag
                    d[f"epsilon_im_{x}{y}"] = eps_im

        return pd.DataFrame(d)

    @staticmethod
    def from_scalar(e, eps_sc, e_units="ev", d=3):
        """Construct a diagonal D-dimensional dielectric function from a
        scalar dielectric function.

        Parameters
        ----------
        eps_sc : complex or array_like
            Scalar dielectric function in units of relative permittivity
            (shape: `(O,)`).
        e : float or array_like, optional
            Energy in `e_units` (array shape: `(O,)`, default: 0.0).
        e_units : str, optional
            Units of `e` (default: "ev").
        d : int
            Dimension of dielectric function (D).
        """

        eps_sc = np.atleast_1d(np.array(eps_sc, dtype=np.complex128))

        if eps_sc.ndim != 1:
            raise ValueError(
                "eps_sc must be a scalar or an array_like with shape (O,)."
            )

        if d < 1 or d > 3:
            raise ValueError("d must be between 1 and 3.")

        eps = np.zeros((eps_sc.shape[0], d, d), dtype=np.complex128)

        for idx in range(d):
            eps[:, idx, idx] = eps_sc

        return DielectricFunction(e, eps, e_units=e_units)

    @staticmethod
    def from_upper_triangle(e, eps_tri, e_units="ev"):
        """Construct a symmetric energy-dependent 3D dielectric
        function from the six independent elements in the upper
        triangle.

        Parameters
        ----------
        e : array_like
            Energy in `e_units` (array shape: `(O,)`).
        eps_sc : array_like
            Independent elements of the dielectric function in units of
            relative permittivity (shape: `(6,)` or `(O, 6)`, order: xx,
            yy, zz, xy, xz, yz).
        e_units : str, optional
            Units of `e` (default: "ev").
        """

        eps_tri = np.atleast_2d(np.array(eps_tri, dtype=np.complex128))

        if eps_tri.ndim != 2 or eps_tri.shape[1] != 6:
            raise ValueError(
                "eps_tri must be an array_like with shape (6,) or (O, 6)"
            )

        eps = np.zeros((eps_tri.shape[0], 3, 3), dtype=np.complex128)

        eps[:, 0, 0] = eps_tri[:, 0]  # xx
        eps[:, 1, 1] = eps_tri[:, 1]  # yy
        eps[:, 2, 2] = eps_tri[:, 2]  # zz
        eps[:, 0, 1] = eps[:, 1, 0] = eps_tri[:, 3]  # xy
        eps[:, 0, 2] = eps[:, 2, 0] = eps_tri[:, 4]  # xz
        eps[:, 1, 2] = eps[:, 2, 1] = eps_tri[:, 5]  # yz

        return DielectricFunction(e, eps, e_units=e_units)
