"""Molecular-orbital contour plotting helpers."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .orbitals import imaginary_fraction, rotate_wavefunction_real


class ContourPlotter:
    """State-free builder for grids of molecular-orbital contour plots."""

    _PLANE_SLICES = {
        # Real-space orbitals follow the FFT array order (z, y, x).
        "xy": (0, "$x$ [a.u.]", "$y$ [a.u.]"),
        "xz": (1, "$x$ [a.u.]", "$z$ [a.u.]"),
        "yz": (2, "$y$ [a.u.]", "$z$ [a.u.]"),
    }

    @staticmethod
    def build_contourplot(
        res,
        system,
        filename=None,
        plane="xy",
        sz=None,
        nrows=1,
        ncols=None,
        levels=9,
        dpi=144,
        ngrid=5,
        tick_rotation=0,
        labels=None,
        plot_energies=True,
        show_imaginary_norm=True,
    ):
        """Generate a grid of contour plots for real-space orbitals.

        Parameters
        ----------
        res : dict
            SCF result containing ``orbc_rs`` and, when energies are plotted,
            ``orbe``.
        system : PeriodicSystem
            Cubic system on whose FFT grid the orbitals are sampled.
        filename : path-like, optional
            Save the completed figure when supplied.
        plane : {'xy', 'xz', 'yz', 'auto'} or sequence of str
            Cartesian plane used for every orbital, or one plane per panel.
            ``'auto'`` selects the central plane with the largest norm for
            each orbital, avoiding accidental plots of nodal planes.
        sz : float, optional
            Half-width of the plot in bohr. Defaults to half the cell edge.
        nrows, ncols : int, optional
            Shape of the subplot grid. ``ncols`` defaults to the number of
            columns needed to show all orbitals.
        levels, dpi, ngrid : int, optional
            Number of contour levels, figure resolution, and axis ticks.
        tick_rotation : float, optional
            Rotation of the x-axis tick labels in degrees.
        labels : sequence of str, optional
            Custom orbital labels.
        plot_energies : bool, optional
            Append orbital energies to panel titles.
        show_imaginary_norm : bool, optional
            Report a residual imaginary norm in panel titles. Disable this to
            hide the diagnostic without changing the plotted orbital.

        Returns
        -------
        tuple
            The Matplotlib ``(figure, axes)`` pair; ``axes`` is always a 2-D
            array.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError("Contour plotting requires matplotlib.") from exc

        orbitals = np.asarray(res.get("orbc_rs"))
        if orbitals.ndim != 4 or orbitals.shape[0] == 0:
            raise ValueError("res['orbc_rs'] must have shape (norb, nz, ny, nx).")
        if len(set(orbitals.shape[1:])) != 1:
            raise ValueError("real-space orbitals must be sampled on a cubic grid.")

        norbitals = len(orbitals)
        if not isinstance(nrows, (int, np.integer)) or nrows < 1:
            raise ValueError("nrows must be a positive integer.")
        if ncols is None:
            ncols = math.ceil(norbitals / nrows)
        if not isinstance(ncols, (int, np.integer)) or ncols < 1:
            raise ValueError("ncols must be a positive integer.")
        if not isinstance(levels, (int, np.integer)) or levels < 3:
            raise ValueError("levels must be an integer of at least three.")
        if not isinstance(ngrid, (int, np.integer)) or ngrid < 2:
            raise ValueError("ngrid must be an integer of at least two.")
        if not np.isfinite(tick_rotation):
            raise ValueError("tick_rotation must be finite.")

        panel_count = min(norbitals, nrows * ncols)
        planes = ContourPlotter._normalise_planes(plane, panel_count)
        if labels is not None and len(labels) < panel_count:
            raise ValueError("labels must contain one value per plotted orbital.")

        energies = res.get("orbe")
        if plot_energies:
            if energies is None or len(energies) < panel_count:
                raise ValueError(
                    "res['orbe'] must contain one value per plotted orbital."
                )
            energies = np.asarray(energies)

        cell_size = float(np.cbrt(system.get_omega()))
        if not np.isfinite(cell_size) or cell_size <= 0:
            raise ValueError("system must have a positive finite cell volume.")
        if sz is None:
            sz = cell_size / 2.0
        if not np.isfinite(sz) or sz <= 0 or sz > cell_size / 2.0:
            raise ValueError("sz must be positive and no larger than half the cell edge.")

        figure, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(2 * ncols + 1, 2 * nrows + 1),
            dpi=dpi,
            squeeze=False,
        )

        for orb_idx, axis in enumerate(axes.flat):
            if orb_idx >= panel_count:
                axis.set_visible(False)
                continue

            orbital = rotate_wavefunction_real(orbitals[orb_idx])
            selected_plane = planes[orb_idx]
            if selected_plane == "auto":
                selected_plane = max(
                    ContourPlotter._PLANE_SLICES,
                    key=lambda value: np.linalg.norm(
                        np.take(
                            orbital.real,
                            orbital.shape[
                                ContourPlotter._PLANE_SLICES[value][0]
                            ] // 2,
                            axis=ContourPlotter._PLANE_SLICES[value][0],
                        )
                    ),
                )
            slice_axis, xlabel, ylabel = ContourPlotter._PLANE_SLICES[
                selected_plane
            ]
            field = np.take(
                orbital.real, orbital.shape[slice_axis] // 2, axis=slice_axis
            )
            field, coordinates = ContourPlotter._prepare_central_plane(
                field, sz, cell_size
            )
            limit = float(np.max(np.abs(field)))
            if not np.isfinite(limit) or limit == 0.0:
                raise ValueError(f"central slice of orbital {orb_idx + 1} is zero.")

            contour_levels = np.linspace(-limit, limit, levels)
            axis.contourf(
                coordinates,
                coordinates,
                field,
                cmap="PiYG",
                levels=contour_levels,
                vmin=-limit,
                vmax=limit,
            )
            axis.contour(
                coordinates,
                coordinates,
                field,
                colors="black",
                levels=contour_levels,
                vmin=-limit,
                vmax=limit,
                linewidths=0.6,
            )

            axis.set_aspect("equal", adjustable="box")
            axis.set_xlim(-sz, sz)
            axis.set_ylim(-sz, sz)
            axis.set_xticks(np.linspace(-sz, sz, ngrid))
            axis.set_yticks(np.linspace(-sz, sz, ngrid))
            axis.tick_params(axis="x", labelrotation=tick_rotation)
            axis.grid(linestyle="--", alpha=0.5)
            axis.set_xlabel(xlabel)
            axis.set_ylabel(ylabel)

            if labels is None:
                title = rf"$\psi_{{{orb_idx + 1}}}$"
            else:
                title = str(labels[orb_idx])
            if plot_energies:
                title += f" ({energies[orb_idx]:.4f} Ht)"
            residual = imaginary_fraction(orbital)
            if show_imaginary_norm and residual > 1e-8:
                title += f"\nIm norm: {residual:.1%}"
            axis.set_title(title)

        figure.tight_layout()
        if filename is not None:
            figure.savefig(Path(filename), dpi=dpi)
        return figure, axes

    @staticmethod
    def _normalise_planes(plane, panel_count):
        if isinstance(plane, str):
            planes = [plane] * panel_count
        else:
            try:
                planes = list(plane)
            except TypeError as exc:
                raise ValueError(
                    "plane must be 'xy', 'xz', 'yz', 'auto', or a sequence."
                ) from exc
            if len(planes) < panel_count:
                raise ValueError("plane must contain one value per plotted orbital.")
            planes = planes[:panel_count]

        valid = {*ContourPlotter._PLANE_SLICES, "auto"}
        invalid = [value for value in planes if value not in valid]
        if invalid:
            raise ValueError(
                "plane values must be 'xy', 'xz', 'yz', or 'auto'."
            )
        return planes

    @staticmethod
    def _prepare_central_plane(field, sz, cell_size):
        """Return a centred field with coordinates at the actual FFT nodes.

        FFT samples cover ``[-L/2, L/2)``.  The positive cell boundary is
        therefore added as a periodic copy of the negative boundary.  Exact
        crop boundaries are interpolated when ``sz`` does not coincide with
        an FFT node.
        """
        npts = field.shape[0]
        fft_coordinates = (
            np.arange(npts, dtype=float) / npts - 0.5
        ) * cell_size
        periodic_coordinates = np.append(fft_coordinates, cell_size / 2.0)
        periodic_field = np.pad(field, ((0, 1), (0, 1)), mode="wrap")

        interior = periodic_coordinates[
            (periodic_coordinates > -sz) & (periodic_coordinates < sz)
        ]
        coordinates = np.concatenate(([-sz], interior, [sz]))
        coordinates = np.unique(coordinates)
        if len(coordinates) < 2:
            raise ValueError("sz is too small for the real-space FFT grid.")

        yy, xx = np.meshgrid(coordinates, coordinates, indexing="ij")
        points = np.column_stack((yy.ravel(), xx.ravel()))
        interpolator = RegularGridInterpolator(
            (periodic_coordinates, periodic_coordinates),
            periodic_field,
            bounds_error=True,
        )
        plotted_field = interpolator(points).reshape(yy.shape)
        return plotted_field, coordinates
