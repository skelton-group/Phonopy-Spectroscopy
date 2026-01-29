# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines for reading and writing files not associated with specific
codes."""


# -------
# Imports
# -------


import numpy as np

from ..structure import Structure

from ..utility.numpy_helper import np_check_shape


# ---------
# XYZ Files
# ---------


def structure_from_xyz(file_path, cubic=True, pad=15.0):
    r"""Read coordinates from an XYZ-format file and return a `Structure`
    object with the molecule placed at the centre of a large unit cell.

    Parameters
    ----------
    file_path : str
        Path to input file.
    cubic : bool, optional
        Selects a cubic (`True`) or rectangular cell (`False`) (default:
        `True`).
    pad : float, optional
        Minimum spacing between atoms in periodic images (default: 15
        Ang).

    Returns
    -------
    struct : Structure
        `Structure` object with the molecule placed at the centre of a
        cubic (`cubic=True`) or rectangular cell (`cubic=False`) with a
        minimum distance of `pad` between atoms in periodic images.
    """

    if pad < 0.0:
        raise ValueError("pad cannot be less than zero.")

    with open(file_path, "r") as f:
        # Read atom count.

        n_at = int(next(f).strip())

        # Skip title line.

        next(f)

        # Read atom data.

        at_syms, at_pos = [], []

        for _ in range(n_at):
            vals = next(f).strip().split()

            at_syms.append(vals[0])
            at_pos.append(np.array([float(v) for v in vals[1:4]]))

        # Determine "bounding box" for molecular structure.

        at_pos = np.array(at_pos, dtype=np.float64)

        p_min = at_pos.min(axis=0)
        p_max = at_pos.max(axis=0)

        m_box = p_max - p_min

        # Determine cell box and lattice parameters.

        c_box = None

        if cubic:
            c_box = (m_box.max() + pad) * np.ones((3,), dtype=np.float64)
        else:
            c_box = m_box + pad

        a, b, c = c_box

        v_latt = np.array(
            [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]], dtype=np.float64
        )

        # Shift atoms to place molecule at the centre of the cell.

        at_pos += ((c_box / 2.0) - p_min)[np.newaxis, :]

        # Build and return a Structure object.

        return Structure(v_latt, at_pos, at_syms, cart_to_frac=True)


# --------------
# XCrysDen files
# --------------


def structure_to_xsf(struct, file_path, comment=None, vecs=None):
    """Write a `Structure` object to an XCrysDen-format file with
    optional vectors attached to the atomic positions.

    Parameters
    ----------
    struct : Structure
        Structure.
    file_path : str
        File path.
    comment : str or None, optional
        Optional comment line to add to the top of the file (default:
        `None`).
    vecs : array_like or None, optional
        Optionally specify a set of vectors (shape: `(3,)`) or `None`
        to attach to each atom.
    """

    if struct.num_atoms == 0:
        raise ValueError("Cannot write an empty structure to an XSF file.")

    if vecs is not None:
        if len(vecs) != struct.num_atoms:
            raise ValueError(
                "If supplied, vecs must have one entry per atom in struct."
            )

        # Better to validate vecs before starting to write the output
        # file.

        vecs_inp = vecs

        vecs = []

        for v in vecs_inp:
            if v is not None:
                v = np.asarray(v, dtype=np.float64)

                if not np_check_shape(v, (3,)):
                    raise ValueError(
                        "Vectors in vecs must be array_like with "
                        "shape (3,)."
                    )

            vecs.append(v)

    with open(file_path, "w") as f:
        if comment is not None:
            # The XSF file in principle allows for multiline comments.

            for line in comment.split("\n"):
                if not line.startswith("#"):
                    line = "# " + line

                f.write(line + "\n")

        f.write("CRYSTAL\n")

        f.write("PRIMVEC\n")

        for v in struct.lattice_vectors:
            f.write("{0: >16.10f}  {1: >16.10f}  {2: >16.10f}\n".format(*v))

        f.write("PRIMCOORD\n")
        f.write("{0} 1\n".format(struct.num_atoms))

        at_num = struct.atomic_numbers()
        at_pos = struct.cartesian_positions()

        fmt_pos = "{0: >3}  {1: >16.10f}  {2: >16.10f}  {3: >16.10f}\n"

        fmt_pos_vec = fmt = (
            "{0: >3}  {1: >16.10f}  {2: >16.10f}  {3: >16.10f} "
            "{4: >14.10f}  {5: >14.10f}  {6: >14.10f}\n"
        )

        if vecs is not None:
            for i, (n, p, v) in enumerate(zip(at_num, at_pos, vecs)):
                if v is not None:
                    f.write(fmt_pos_vec.format(n, *p, *v))
                else:
                    f.write(fmt_pos.format(n, *p))

        else:
            for n, p in zip(at_num, at_pos):
                f.write(fmt_pos.format(n, *p))
