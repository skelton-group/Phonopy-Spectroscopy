# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Class for storing and interpolating energy-dependent Raman tensors."""

# -------
# Imports
# -------

import numpy as np

from scipy.interpolate import interp1d

from ..constants import ZERO_TOLERANCE

from ..utility.numpy_helper import (
    np_asarray_copy,
    np_readonly_view,
    np_check_shape,
    np_expand_dims,
    np_discard_imag_if_real,
)

# ------------------
# RamanTensors class
# ------------------


class RamanTensors:
    """Store and interpolate energy-dependent Raman tensors."""

    def __init__(self, r_t, e=0.0):
        """Create a new `RamanTensors` instance.

        Parameters
        ----------
            r_t : array_likenp_asarray_copy
                Set of `N` Raman tensors (shape: `(N, 3, 3)`)
                or `N` Raman tensors at `M` photon energies
                (shape: `(N, M, 3, 3)`).
            e : array_like, optional
                For energy-dependent Raman tensors, `M` photon energies
                in eV. `e` must be monotinically increasing with
                `e[0] == 0.0`. When `e` is not specified, the `N` Raman
                tensors are assumed to be calculated at E = 0.
        """

        r_t, _ = np_expand_dims(
            np_asarray_copy(r_t, dtype=np.complex128),
            (None, None, 3, 3),
            expand_order=(1, 0),
        )

        if np.ndim(e) == 0:
            e = np.array([e], dtype=np.float64)

        e = np_asarray_copy(e, dtype=np.float64)

        _, n_e, _, _ = r_t.shape

        if not np_check_shape(e, (n_e,)):
            raise ValueError("e must have shape (M,).")

        # Negative photon energies are unphysical.

        if (e < 0.0).any():
            raise ValueError("Photon energies cannot be negative.")

        # For energy-dependent calculations, ensure the E are in
        # increasing order (needed for interpolation). This, plus the
        # condition above, also ensures that all E are positive.

        if len(e) > 1 and (e[1:] <= e[:-1]).any():
            raise ValueError("e must be monotonically increasing.")

        # Using complex instead of real numbers in calculations can
        # incur a performance penalty. If the Raman tensors are real,
        # drop the imaginary part and convert them to np.float64.

        r_t = np_discard_imag_if_real(r_t)

        self._e = e
        self._r_t = r_t

    @property
    def energies(self):
        """numpy.ndarray : Photon energies in eV (shape: `(M,)`)."""
        return np_readonly_view(self._e)

    @property
    def raman_tensors(self):
        """numpy.ndarray : Raman tensors in Ang^2 / sqrt(amu) (shape:
        `(N, M, 3, 3)`)."""
        return np_readonly_view(self._r_t)

    @property
    def num_modes(self):
        """int : Number of Raman tensors (`N`)."""
        return self._r_t.shape[0]

    @property
    def num_energies(self):
        "int : Number of photon energies (`M`)."
        return self._r_t.shape[1]

    @property
    def is_energy_dependent(self):
        """bool : `True` if the Raman tensors are energy-dependent,
        otherwise `False`."""
        return len(self._e) > 1

    def get_tensors_at_energy(self, e):
        """Get Raman tensors at a specific photon energy.

        Parameters
        ----------
        e : float
            Photon energy in eV.

        Returns
        -------
        r_t : numpy.ndarray
            Raman tensors at supplied energy (shape: `(N, 3, 3)`).
        """

        # If the requested E matches any of the calculated E, return
        # that set of tensors.

        for idx, e_calc in enumerate(self._e):
            if np.abs(e - e_calc) <= ZERO_TOLERANCE:
                return self._r_t[:, idx, :, :].reshape((-1, 3, 3))

        # If not, and if E falls within the calculated range,
        # interpolate.

        if e < self._e.min() or e > self._e.max():
            raise RuntimeError(
                "Requested E = {0:.3f} is outside the calculated "
                "E = {1:.3f} to {2:.3f}.".format(
                    e, self._e.min(), self._e.max()
                )
            )

        r_t_interp = np.zeros(
            (self._r_t.shape[0], 3, 3), dtype=self._r_t.dtype
        )

        for i, t in enumerate(self._r_t):
            for j in range(3):
                for k in range(3):
                    f = interp1d(self._e, t[:, j, k], kind="linear")
                    r_t_interp[i, j, k] = f(e)

        return r_t_interp

    def to_dict(self):
        """Return the internal data as a dictionary of native Python
        types for serialisation.

        Returns
        -------
        d : dict
            Dictionary structure containing internal data as native
            Python types.
        """

        # JSON cannot store complex numbers, so if we have complex Raman
        # tensors we need to serialise the real and imaginary parts
        # separately.

        r_t_dict = {"real": self._r_t.real.tolist()}

        if np.iscomplexobj(self._r_t):
            r_t_dict["imag"] = self._r_t.imag.tolist()

        return {"energies": self._e.tolist(), "raman_tensors": r_t_dict}

    @staticmethod
    def from_dict(d):
        """Create a new `RamanTensors` instance from a dictionary
        generated by `RamanTensors.to_dict()`.

        Parameters
        ----------
        d : dict
            Dictionary generated by `to_dict()`.

        Returns
        -------
        r_t : RamanTensors
            `RamanTensors` object constructed from the data in `d`.
        """

        # The "raman_tensors" key may contain only a "real" key, or may
        # contain both "real" and "imag" keys.

        r_t_real = np.asarray(d["raman_tensors"]["real"], dtype=np.float64)

        r_t_imag = None

        if "imag" in d["raman_tensors"]:
            r_t_imag = np.asarray(d["raman_tensors"]["imag"], dtype=np.float64)

            if not np.equal(np.shape(r_t_imag), np.shape(r_t_real)).all():
                raise ValueError(
                    '"raman_tensors" key contains "real" and "imag" '
                    "keys with different array shapes."
                )

        r_t = None

        if r_t_imag is None:
            r_t = np_asarray_copy(r_t_real, dtype=np.float64)
        else:
            r_t = np.zeros_like(r_t_real, dtype=np.complex128)

            r_t.real = r_t_real
            r_t.imag = r_t_imag

        return RamanTensors(r_t, e=d["energies"])
