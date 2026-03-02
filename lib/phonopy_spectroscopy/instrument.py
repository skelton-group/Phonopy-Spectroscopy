# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Classes for defining instrument geometries and polarisations."""


# -------
# Imports
# -------

from itertools import product

import numpy as np

from .constants import ZERO_TOLERANCE

from .utility.geometry import (
    parse_direction,
    rotation_matrix_from_axis_angle,
    rotation_matrix_from_vectors,
)

from .utility.numpy_helper import (
    np_asarray_copy,
    np_readonly_view,
    np_check_shape,
    np_expand_dims,
)

from .utility.quadrature import unit_circle_quad_rule


# --------------
# Geometry class
# --------------


class Geometry:
    """Represent a measurement geometry in a polarised Raman
    experiment."""

    def __init__(self, i_dir, c_dir):
        """Create a new instance of the `Geometry` class.

        Parameters
        ----------
        i_dir, c_dir : array_like or str
            Incident and collection directions specified as Cartesian
            directions or vectors.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `i_dir`
            and `c_dir`.
        """

        # parse_direction will validate the direction specifiers and
        # normalise the vectors.

        self._i_dir = parse_direction(i_dir, norm=True)
        self._c_dir = parse_direction(c_dir, norm=True)

    @property
    def incident_direction(self):
        """numpy.ndarray : Incident light direction (shape: `(3,)`)."""
        return np_readonly_view(self._i_dir)

    @property
    def collection_direction(self):
        """numpy.ndarray : Collection direction (shape: `(3,)`)."""
        return np_readonly_view(self._c_dir)

    def check_incident_polarisations(self, i_pols):
        """Check whether incident light polarisations are valid for the
        geometry.

        Parameters
        ----------
        i_pols : Polarisation or array_like
            Polarisation(s) to check.

        Returns
        -------
        valid : bool
            `True` if the polarisation(s) is/are perpendicular to the
            incident light direction, otherwise `False`.
        """

        i_pols, _ = np_expand_dims(np.asarray(i_pols, dtype=object), (None,))

        for i_pol in i_pols:
            if not i_pol.check_perpendicular(self._i_dir):
                return False

        return True

    def check_scattered_polarisations(self, s_pols):
        """Check whether scattered light polarisations are valid for the
        geometry.

        Parameters
        ----------
        s_pols : Polarisation or array_like
            Polarisation(s) to check.

        Returns
        -------
        valid : bool
            `True` if polarisation(s) is/are perpendicular to the
            collection axis, otherwise `False`.
        """

        s_pols, _ = np_expand_dims(np.asarray(s_pols, dtype=object), (None,))

        for s_pol in s_pols:
            if not s_pol.check_perpendicular(self._c_dir):
                return False

        return True


# ------------------
# Polarisation class
# ------------------


class Polarisation:
    """Represent a polarisation in a polarised Raman measurement."""

    def __init__(self, v, w=None):
        """Create a new instance of the `Polarisation` class.

        Parameters
        ----------
        v : array_like
            3D polarisation vector (shape: `(3,)`) or vectors (shape:
            `(N, 3)`).
        w : array_like or None, optional
            Weights for summing/averaging multiple `v` (shape: `(N,)`).

        See Also
        --------
        Polarisation.from_direction :
            Define a polarisation from a direction.
        Polarisation.from_angles :
            Define (a) polarisation(s) as (a) rotation(s) about an axis.
        Polarsation.from_rotation :
            Define polarisations for an angle rotation about an axis.
        Polarisation.integration :
            Define an integration over polarisations.
        Polarisation.cross_to :
            Define polarisation(s) cross to (an)other(s).
        Polarisation.sum_parallel_cross_to :
            Define (a) sum(s) of polarisation(s) parallel and cross to
            (an)other.

        Notes
        -----
        For most use cases, it is likely more convenient to create
        `Polarisation` objects using the static methods on this class
        than to instantiate them directly (see above).
        """

        v, _ = np_expand_dims(
            np_asarray_copy(v, dtype=np.complex128), (None, 3)
        )

        # For most "routine" calculations the polarisation vectors are
        # real - if so, drop the imaginary part and convert to
        # np.float64 for performance.

        if not np.iscomplex(v).any():
            v = np.array(v.real, dtype=np.float64)

        if w is not None:
            w = np_asarray_copy(w, dtype=np.float64)

            if len(w) != v.shape[0]:
                raise ValueError("w must be an array_like with shape (N,).")
        else:
            if v.shape[0] > 1:
                raise ValueError(
                    "w must be specified for multiple polarisation vectors."
                )

            w = np.array([1.0], dtype=np.float64)

        # Check vectors are non-zero and normalised.

        if not np.allclose(
            np.linalg.norm(v, axis=1), 1.0, atol=ZERO_TOLERANCE
        ):
            raise ValueError(
                "Polarisation vectors must be non-zero and normalised."
            )

        self._v = v
        self._w = w

    @property
    def vectors(self):
        """numpy.ndarray : Polarisation vectors (shape: `(N, 3)`)."""
        return np_readonly_view(self._v)

    @property
    def weights(self):
        """numpy.ndarray : Polarisation vector weights (shape: `(N,)`)."""
        return np_readonly_view(self._w)

    @property
    def num_vecs_w(self):
        """int : Number of weights/vectors."""
        return self._v.shape[0]

    @property
    def is_complex(self):
        """bool : `True` if any of `vectors` are complex, otherwise
        `False`."""
        return np.iscomplex(self._v).any()

    def iter_v_w(self):
        """Iterate over polarisation vectors and weights.

        Yields
        ------
        v_w : tuple of (numpy.ndarray, float)
            Polarisation vector and weight.
        """

        for v_w in zip(self._v, self._w):
            yield v_w

    def check_perpendicular(self, axis):
        """Check polarisation vectors are perpendicular to an axis.

        Parameters
        ----------
        axis : array_like or str
            Axis to check.

        Returns
        -------
        perp : bool
            `True` if all polarisation vectors are perpendicular to
            axis, otherwise `False`.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `axis`.
        """

        axis = parse_direction(axis)

        for v in self._v:
            # cos(\theta) = 0 for perpendicular vectors, so it does not
            # matter whether axis is normalised.

            if np.abs(np.dot(axis, v)) > ZERO_TOLERANCE:
                return False

        return True

    def combine_with(self, other):
        """Return the product of the polarisation weights and vectors in
        in this `Polarisation` object with those in another.

        Parameters
        ----------
        other : Polarisation
            `Polarisation` object to combine with.

        Returns
        -------
        prod : list of tuples of (numpy.ndarray, numpy.ndarray, float)
            List of `(v_s, v_o, w)` tuples containing the vectors from
            this ("self") and the other polarisations (shape: `(3,)`),
            and the combined weight `w`.
        """

        return list(self.combine_with_iter(other))

    def combine_with_iter(self, other):
        """Iterate over the product of the polarisation weights and
        vectors in in this `Polarisation` object with those in another.

        Parameters
        ----------
        other : Polarisation
            `Polarisation` object to combine with.

        Yields
        ------
        y : tuple of (numpy.ndarray, numpy.ndarray, float)
            Tuples of `(v_s, v_o, w)` containing the vectors from this
            ("self") and the other polarisations (shape: `(3,)`), and
            the combined weight `w`.
        """

        for (v_s, w_s), (v_o, w_o) in product(
            self.iter_v_w(), other.iter_v_w()
        ):
            yield (v_s, v_o, w_s * w_o)

    @staticmethod
    def from_direction(dirn):
        """Define a polarisation from a direction.

        Parameters
        ----------
        v : array_like or str
            Direction.

        Returns
        -------
        p : Polarisation
            A `Polarisation` object for direction `dirn`.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `dirn`.
        """

        return Polarisation(parse_direction(dirn))

    @staticmethod
    def from_angles(axis, angles):
        """Define (a) polarisation(s) as (a) rotation(s) about an axis.

        Parameters
        ----------
        axis : array_like or str
            Axis of rotation.
        angles : float or array_like
            Angle(s) of rotation.

        Returns
        -------
        pols : Polarisation or numpy.ndarray
            `Polarisation` or array of `Polarisation` (same shape as
            `angles`).

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `axis`.
        """

        axis = parse_direction(axis)

        angles, n_dim_add = np_expand_dims(
            np.asarray(angles, dtype=np.float64), (None,)
        )

        # Obtain an initial polarisation vector v_0 perpendicular to the
        # rotation axis by determining the rotation matrix that rotates
        # the axis to +z and applying it to +x.

        r = rotation_matrix_from_vectors(axis, parse_direction("z"))
        v_0 = np.matmul(r, parse_direction("x"))

        # Generate a sequence of polarisation vectors by rotating v_0
        # around axis.

        vecs = [
            np.matmul(rotation_matrix_from_axis_angle(axis, theta), v_0)
            for theta in angles
        ]

        pols = np.asarray(
            [Polarisation([v], [1.0]) for v in vecs], dtype=object
        )

        return pols if n_dim_add == 0 else pols[0]

    @staticmethod
    def from_rotation(axis, start=0.0, end=360.0, step=2.5):
        """Define polarisations for an angle rotation about an axis.

        Parameters
        ----------
        axis : array_like or str
            Axis of rotation.
        start, end, step : float, optional
            Start/end angle and angle step in degrees (defaults:
            start = 0.0, end = 360.0, step = 2.5).

        Returns
        -------
        pols : numpy.ndarray
            `Polarisation` objects for each step in the angle rotation.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `axis`.
        """

        angles = np.arange(start, end + step / 10.0, step)

        if len(angles) == 0:
            raise ValueError(
                "No angles between start = {0:.2f} -> end = {1:.2f} "
                "with step = {2:.2f}.".format(start, end, step)
            )

        return Polarisation.from_angles(axis, angles)

    @staticmethod
    def integration(axis, n=16):
        """Define an integration over polarisations perpendicular to an
        axis using a circle quadrature rule.

        Parameters
        ----------
        axis : array_like or str
            Axis to which polarisation vectors should be perpendicular.
        n : int
            Number of points for the circle quadrature rule.

        Returns
        -------
        p : Polarisation
            A `Polarisation` object for the integration.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `axis`.
        """

        axis = parse_direction(axis)

        # Obtain a rotation matrix for rotating vectors in the (x, y)
        # plane to be perpendicular to axis by determining the rotation
        # matrix that rotates the axis to +z.

        r = rotation_matrix_from_vectors(axis, parse_direction("z"))

        # Get vectors and weights for circle rule.

        vecs, w = unit_circle_quad_rule(n, ret="vectors")

        # Double transpose to convert vecs to/from column format.

        return Polarisation(np.matmul(r, vecs.T).T, w)

    @staticmethod
    def cross_to(pol, axis, rot_dir=1.0):
        """Define polarisation(s) cross to (an)other(s).

        Parameters
        ----------
        pol : Polarisation or array_like
            Polarisation(s) to cross.
        axis : array_like or str
            Axis to cross polarisation.
        rot_dir : int
            Sign of rotation (+ve = anticlockwise, -ve = clockwise;
            default: +1.0).

        Returns
        -------
        pols : Polarisation or numpy.ndarray
            `Polarisation` or array of `Polarisation` (same shape as
            `pol`).

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `axis`.
        """

        axis = parse_direction(axis)

        pol, n_dim_add = np_expand_dims(np.asarray(pol, dtype=object), (None,))

        for p in pol:
            if not p.check_perpendicular(axis):
                raise ValueError("axis must be perpendicular to (all) pol.")

        r = rotation_matrix_from_axis_angle(axis, np.sign(rot_dir) * 90.0)

        pol_cross = np.zeros_like(pol)

        for i, p in enumerate(pol):
            pol_cross[i] = Polarisation(np.matmul(r, p.vectors.T).T, p.weights)

        return pol_cross if n_dim_add == 0 else pol_cross[0]

    @staticmethod
    def sum_parallel_cross_to(pol, axis, rot_dir=1.0):
        """Define (a) sum(s) of polarisation(s) parallel and cross to
        (an)other.

        Parameters
        ----------
        pol : Polarisation or array_like
            Polarisation(s) to sum parallel/cross.
        axis : array_like or str
            Axis to cross polarisation.
        rot_dir : int
            Sign of rotation for cross polarisation (+ve =
            anticlockwise, -ve = clockwise; default: +1.0).

        Returns
        -------
        pols : Polarisation or numpy.ndarray
            `Polarisation` or array of `Polarisation (same shape as
            `pol`).

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `axis`.
        """

        pol, n_dim_add = np_expand_dims(np.asarray(pol, dtype=object), (None,))

        pol_cross = Polarisation.cross_to(pol, axis, rot_dir=rot_dir)

        p_sum = np.zeros_like(pol)

        for i, (p_par, p_cross) in enumerate(zip(pol, pol_cross)):
            vecs, w = [], []

            for (v_1, w_1), (v_2, w_2) in zip(
                p_par.iter_v_w(), p_cross.iter_v_w()
            ):
                vecs.append(v_1)
                vecs.append(v_2)

                w.append(w_1)
                w.append(w_2)

            p_sum[i] = Polarisation(vecs, w)

        return p_sum if n_dim_add == 0 else p_sum[0]
