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

# --------------
# Geometry class
# --------------


class Geometry:
    """Represent a measurement geometry."""

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

        # parse_direction will validate the direction specifiers,
        # normalise the vectors. This also ensures it is copied.

        i_dir = parse_direction(i_dir)
        c_dir = parse_direction(c_dir)

        i_dir.flags.writeable = False
        c_dir.flags.writeable = False

        self._i_dir = i_dir
        self._c_dir = c_dir

    @property
    def incident_direction(self):
        """numpy.ndarray : Incident light direction (shape: `(3,)`)."""
        return self._i_dir

    @property
    def collection_direction(self):
        """numpy.ndarray : Collection direction (shape: `(3,)`)."""
        return self._c_dir

    @property
    def raman_nac_direction(self):
        """numpy.ndarray : Approach direction q in Raman measurements
        for computing non-analytical corrections to the dynamical
        matrix."""

        q = self._i_dir - self._c_dir

        if (np.abs(q) > ZERO_TOLERANCE).any():
            return q / np.linalg.norm(q)
        else:
            return q

    @staticmethod
    def conventional_backscattering():
        """Geometry for a conventional backscattering measurement."""

        return Geometry("+z", "-z")

    @staticmethod
    def conventional_off_axis_transmission(theta, n_f=1.0, n_b=1.0):
        r"""Geometry for a conventional off-axis transmission
        measurement, optionally accounting for the refractive indices
        of the "front" (incident) and "back" (exit) media.

        Parameters
        ----------
        theta : float
            Incidence angle in degrees.
        n_f, n_b : float, optional
            Refractive indices of the front and back media (default: 1.0
            = vacuum ~ air).

        Returns
        -------
        geom : Geometry
            Measurement geometry.

        Notes
        -----
        This assumes a conventional geometry where light is incident
        in the xz plane at an angle \theta to the z axis.
        """

        t = np.radians(theta)

        q_x = n_f * np.sin(t)
        q_i_z = n_f * np.cos(t)

        # At oblique angles the z component of the transmitted
        # wavevector may be complex.

        q_t_z = np.sqrt(n_b**2 - q_x**2, dtype=np.complex128)

        return Geometry([q_x, 0.0, q_i_z], [q_x, 0.0, q_t_z])

    @staticmethod
    def conventional_off_axis_reflectivity(theta, n_f=1.0):
        r"""Geometry for a conventional off-axis reflectivity
        measurement, optionally accounting for the refractive index of
        the "front" (indicent) medium.

        Parameters
        ----------
        theta : float
            Incidence angle in degrees.
        n_f : float, optional
            Refractive indices of the front medium (default: 1.0
            = vacuum ~ air).

        Returns
        -------
        geom : Geometry
            Measurement geometry.

        Notes
        -----
        This assumes a conventional geometry where light is incident
        in the xz plane at an angle \theta to the z axis.
        """

        t = np.radians(theta)

        q_x = n_f * np.sin(t)
        q_z = n_f * np.cos(t)

        return Geometry([q_x, 0.0, q_z], [q_x, 0.0, -1.0 * q_z])


# ----------------------
# PolarisationBase class
# ----------------------


class PolarisationBase:
    """Represent a D-dimensional polarisation or weighted sum/average of
    polarisations."""

    def __init__(self, v, w=None):
        """Create a new instance of the `PolarisationBase` class.

        Parameters
        ----------
        v : array_like
            D-dimensional polarisation vectors (shape: `(D,)` or
            `(N, D)`).
        w : array_like or None, optional
            Weights for summing/averaging `v` (shape: `(N,)`; default:
            equal weight).
        """

        v = np.atleast_2d(np.array(v, dtype=np.complex128, copy=True))

        if v.ndim != 2:
            raise ValueError(
                "v must be an array_like with shape (D,) or (N, D)."
            )

        # For most "routine" calculations the polarisation vectors are
        # real, and discarding the imaginary part improves performance.

        v = np.real_if_close(v, tol=ZERO_TOLERANCE)

        # Check vectors are non-zero and normalised.

        if not np.allclose(
            np.linalg.norm(v, axis=1), 1.0, atol=ZERO_TOLERANCE
        ):
            raise ValueError("All v must be non-zero and normalised.")

        if w is not None:
            w = np.atleast_1d(np.array(w, dtype=np.float64, copy=True))

            if w.ndim != 1 or len(w) != len(v):
                raise ValueError("w must be an array_like with shape (N,).")
        else:
            w = np.ones((len(v),), dtype=np.float64)

        v.flags.writeable = False
        w.flags.writeable = False

        self._v = v
        self._w = w

    @property
    def vectors(self):
        """numpy.ndarray : Polarisation vectors (shape: `(N, D)`)."""
        return self._v

    @property
    def weights(self):
        """numpy.ndarray : Polarisation vector weights (shape: `(N,)`)."""
        return self._w

    @property
    def dimension(self):
        """int : Dimension of polarisation vectors."""
        return self._v.shape[1]

    def iter_v_w(self):
        """Iterate over polarisation vectors and weights.

        Yields
        ------
        v_w : tuple of (numpy.ndarray, float)
            Polarisation vector and weight.
        """

        for v_w in zip(self._v, self._w):
            yield v_w

    def combine_with(self, other):
        """Combine with another `Polarisation` object to generate pairs
        of polarisation vectors and combined weights.

        Parameters
        ----------
        other : Polarisation
            `Polarisation` object to combine with.

        Returns
        -------
        prod : list of tuples of (numpy.ndarray, numpy.ndarray, float)
            List of `(v_s, v_o, w)` tuples containing pairs of
            polarisation vectors (shape: `(3,)`) and weights.
        """

        return list(self.combine_with_iter(other))

    def combine_with_iter(self, other):
        """Combine with another `Polarisation` object and iterate over
        pairs of polarisation vectors and combined weights.

        Parameters
        ----------
        other : Polarisation
            `Polarisation` object to combine with.

        Yields
        ------
        y : tuple of (numpy.ndarray, numpy.ndarray, float)
            `(v_s, v_o, w)` tuples containing pairs of
            polarisation vectors (shape: `(3,)`) and weights.
        """

        if other.dimension != self._v.shape[1]:
            raise ValueError(
                "Polarisations can only be combined if they have the "
                "same dimension."
            )

        for (v_s, w_s), (v_o, w_o) in product(
            self.iter_v_w(), other.iter_v_w()
        ):
            yield (v_s, v_o, w_s * w_o)


# ------------------
# Polarisation class
# ------------------


class Polarisation(PolarisationBase):
    """Represent a 3D polarisation or weighted sum/average of
    polarisations."""

    def __init__(self, v, w=None):
        """Create a new instance of the `Polarisation` class.

        Parameters
        ----------
        v : array_like
            3D polarisation vectors (shapes: `(3,)` or `(N, 3)`).
        w : array_like or None, optional
            Weights for summing/averaging `v` (shape: `(N,)`; default:
            equal weight).

        See Also
        --------
        Polarisation.from_direction :
            Define a polarisation from a direction.
        Polarisation.from_angles :
            Define polarisations as rotations about an axis.
        Polarsation.from_rotation :
            Define polarisations from an angle rotation about an axis.
        Polarisation.cross_to :
            Define a cross polarisation.
        Polarisation.sum_parallel_cross_to :
            Define a sum of parallel and cross polarisations.

        Notes
        -----
        For most use cases, it is likely more convenient to create
        `Polarisation` objects using the static methods on this class
        than to instantiate them directly.
        """

        super(Polarisation, self).__init__(v, w=w)

        if self._v.shape[1] != 3:
            raise ValueError(
                "Polarisation objects must be initialised with 3D vectors."
            )

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
            `axis`, otherwise `False`.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `axis`.
        """

        axis = parse_direction(axis)

        for v in self._v:
            # cos(\theta) = 0 for perpendicular vectors

            if np.abs(np.dot(axis, v)) > ZERO_TOLERANCE:
                return False

        return True

    @staticmethod
    def conventional_horizontal(i_dir="z"):
        r"""Conventional "horizontal" polarisation.

        Parameters
        ----------
        i_dir : array_like or str
            Incident light direction.

        Returns
        -------
        p : Polarisation
            `Polarisation` object equivalent to "+x" in a conventional
            transmission or backscattering geometry.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `i_dir`.

        Notes
        -----
        Horizontal polarisation :math:`\hat{x}^\prime` is defined as
        :math:`\hat{x}` rotated to be perpendicular to `i_dir`.
        """

        # Rotation matrix that rotates z onto i_dir.

        r = rotation_matrix_from_vectors("z", i_dir)

        return Polarisation.from_direction(r @ parse_direction("x"))

    @staticmethod
    def unpolarised_incident(i_dir="z"):
        r"""Unpolarised incident light.

        Parameters
        ----------
        i_dir : array_like or str
            Incident light direction.

        Returns
        -------
        p : Polarisation
            `Polarisation` equivalent to :math:`(x + y) / 2` in a
            conventional transmission or backscattering geometry.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `i_dir`.

        Notes
        -----
        Unpolarised incident light is defined as:

        .. math::

            \frac{1}{2} \hat{x}^\prime + \frac{1}{2} \hat{y}^\prime

        where :math:`\hat{x}^\prime` and :math:`\hat{y}^\prime` are the
        Cartesian vectors :math:`\hat{x}` and :math:`\hat{y}` rotated
        to be perpendicular to `i_dir`.
        """

        r = rotation_matrix_from_vectors("z", i_dir)

        x_p = np.matmul(r, parse_direction("x"))
        y_p = np.matmul(r, parse_direction("y"))

        return Polarisation([x_p, y_p], [0.5, 0.5])

    @staticmethod
    def unpolarised_collected(c_dir="z"):
        r"""Unpolarised collected light.

        Parameters
        ----------
        c_dir : array_like or str
            Collected light direction.

        Returns
        -------
        p : Polarisation
            `Polarisation` equivalent to :math:`x + y` in a conventional
            transmission geometry.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `c_dir`.

        Notes
        -----
        Unpolarised collected light is defined as:

        .. math::

            \hat{x}^\prime + \hat{y}^\prime

        where :math:`\hat{x}^\prime` and :math:`\hat{y}^\prime` are the
        Cartesian vectors :math:`\hat{x}` and :math:`\hat{y}` rotated
        to be perpendicular to `c_dir`.
        """

        r = rotation_matrix_from_vectors("z", c_dir)

        x_p = np.matmul(r, parse_direction("x"))
        y_p = np.matmul(r, parse_direction("y"))

        return Polarisation([x_p, y_p], [1.0, 1.0])

    @staticmethod
    def from_direction(dirn):
        """Polarisation from a direction.

        Parameters
        ----------
        v : array_like or str
            Direction.

        Returns
        -------
        p : Polarisation
            `Polarisation` for `dirn`.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `dirn`.
        """

        return Polarisation(parse_direction(dirn))

    @staticmethod
    def from_angles(axis, angles):
        """Polarisations from rotations about an axis.

        Parameters
        ----------
        axis : array_like or str
            Rotation axis.
        angles : scalar or array_like
            Rotation angle(s) in degrees.

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

        angles = np.asanyarray(angles, dtype=np.float64)

        ret_scalar = angles.ndim == 0

        angles = np.radians(np.atleast_1d(angles))

        # Obtain an initial polarisation vector v_0 perpendicular to the
        # rotation axis by determining the rotation matrix that rotates
        # the axis to +z and applying it to +x.

        r = rotation_matrix_from_vectors(axis, parse_direction("z"))
        v_0 = r @ parse_direction("x")

        # Generate a sequence of polarisation vectors by rotating v_0
        # around axis.

        pols = np.zeros_like(angles, dtype=object)

        for idx in np.ndindex(*angles.shape):
            r = rotation_matrix_from_axis_angle(axis, angles[idx])
            pols[idx] = Polarisation(r @ v_0)

        return pols[0] if ret_scalar else pols

    @staticmethod
    def from_rotation(axis, start=0.0, end=360.0, step=2.5):
        """Polarisations from a rotation about an axis.

        Parameters
        ----------
        axis : array_like or str
            Rotation axis.
        start, end, step : float, optional
            Start/end angle and angle step in degrees (defaults:
            start = 0.0, end = 360.0, step = 2.5).

        Returns
        -------
        pols : numpy.ndarray
            `Polarisation`s for each angle.

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
    def cross_to(pol, axis, rot_dir=1.0):
        """Perpendicular (cross) polarisation.

        Parameters
        ----------
        pol : Polarisation or array_like
            Polarisations to cross.
        axis : array_like or str
            Axis to cross polarisations.
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

        pol = np.asanyarray(pol, dtype=object)

        ret_scalar = pol.ndim == 0

        pol = np.atleast_1d(pol)

        for p in pol.flat:
            if not p.check_perpendicular(axis):
                raise ValueError("axis must be perpendicular to (all) pol.")

        r = rotation_matrix_from_axis_angle(
            axis, np.sign(rot_dir) * np.radians(90.0)
        )

        pol_cross = np.zeros_like(pol)

        for idx in np.ndindex(*pol.shape):
            p = pol[idx]
            pol_cross[idx] = Polarisation(p.vectors @ r.T, p.weights)

        return pol_cross[0] if ret_scalar else pol_cross

    @staticmethod
    def sum_parallel_cross_to(pol, axis, rot_dir=1.0):
        """Sum of parallel and perpendicular (cross) polarisation.

        Parameters
        ----------
        pol : Polarisation or array_like
            Polarisations to sum parallel/cross.
        axis : array_like or str
            Axis to cross polarisations.
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

        pol = np.asanyarray(pol, dtype=object)

        ret_scalar = pol.ndim == 0

        pol = np.atleast_1d(pol)

        # Generate cross polarisations.

        pol_cross = Polarisation.cross_to(pol, axis, rot_dir=rot_dir)

        p_sum = np.zeros_like(pol, dtype=object)

        for idx in np.ndindex(*pol.shape):
            p_par = pol[idx]
            p_per = pol_cross[idx]

            p_sum[idx] = Polarisation(
                np.concatenate((p_par.vectors, p_per.vectors), axis=0),
                np.concatenate((p_par.weights, p_per.weights), axis=0),
            )

        return p_sum[0] if ret_scalar else p_sum


# --------------------
# SPPolarisation class
# --------------------


class SPPolarisation(PolarisationBase):
    """Represent a polarisation or weighted sum/average of polarisations
    defined relative to a plane of incidence."""

    def __init__(self, v, w=None):
        """Create a new instance of the `SPPolarisation` class.

        Parameters
        ----------
        v : array_like
            2D polarisation vectors (shapes: `(2,)` or `(N, 2)`).
        w : array_like or None, optional
            Weights for summing/averaging `v` (shape: `(N,)`; default:
            equal weight).
        """

        super(SPPolarisation, self).__init__(v, w=w)

        if self._v.shape[1] != 2:
            raise ValueError(
                "SPPolarisation objects must be initialised with 2D "
                "vectors."
            )

    def to_polarisation(self, q, n="+z"):
        r"""Convert the s/p polarisation to a standard (3D) polarisation
        for a given wavevector `q` and surface normal `n`.

        Parameters
        ----------
        q : array_like or str
            Incident wavevector.
        n : array_like or str, optional
            Normal vector (default: "+z").

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `q` and
        `n`.

        Notes
        -----
        If `q` and `n` are collinear, the conventional s and p
        directions :math:`\hat{s} = \hat{y}` and
        :math:`\hat{p} = \hat{x}` are assumed. Otherwise, the s and p
        directions are given by:

        .. math::

            \hat{s} = \hat{n} \times \hat{q}

        .. math::

            \hat{p} = \hat{s} \times \hat{q}
        """

        q = parse_direction(q)
        n = parse_direction(n)

        s_dirn, p_dirn = None, None

        if (np.abs(np.dot(q, n)) - 1.0) < ZERO_TOLERANCE:
            # q and n are parallel -> use default convention.

            s_dirn = parse_direction("y")
            p_dirn = parse_direction("x")
        else:
            s_dirn = np.cross(n, q)
            p_dirn = np.cross(s_dirn, q)

        v_3d = [v[0] * s_dirn + v[1] * p_dirn for v in self._v]

        return Polarisation(v_3d, w=self._w)

    @staticmethod
    def unpolarised_incident():
        """Unpolarised incident light `(s + p) / 2`.

        Returns
        -------
        pol : SPPolarisation
            Polarisation.
        """

        return SPPolarisation([[1.0, 0.0], [0.0, 1.0]], w=[0.5, 0.5])

    @staticmethod
    def unpolarised_collected():
        """Unpolarised collected light `s + p`.

        Returns
        -------
        pol : SPPolarisation
            Polarisation.
        """

        return SPPolarisation([[1.0, 0.0], [0.0, 1.0]], w=[1.0, 1.0])
