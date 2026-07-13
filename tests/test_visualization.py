import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pypwdft.psystem import PeriodicSystem
from pypwdft.visualization import ContourPlotter


@pytest.fixture
def contour_data():
    system = PeriodicSystem(8.0, ecut=2.0)
    coordinates = np.linspace(-1.0, 1.0, 8, endpoint=False)
    z, y, x = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    orbitals = np.asarray([
        x + 0.2 * y,
        (y + 0.2 * z) * np.exp(0.4j),
        z + 0.2 * x,
    ])
    result = {"orbc_rs": orbitals, "orbe": [-0.5, -0.2, 0.1]}
    return system, result


def test_contour_plotter_builds_plane_specific_grid(contour_data, tmp_path):
    system, result = contour_data
    output = tmp_path / "orbitals.png"

    figure, axes = ContourPlotter.build_contourplot(
        result,
        system,
        output,
        plane=["xy", "xz", "yz"],
        nrows=2,
        ncols=2,
        labels=["x", "y", "z"],
    )

    assert axes.shape == (2, 2)
    assert not axes[1, 1].get_visible()
    assert axes[0, 0].get_xlabel() == "$x$ [a.u.]"
    assert axes[0, 1].get_ylabel() == "$z$ [a.u.]"
    assert axes[1, 0].get_xlabel() == "$y$ [a.u.]"
    assert axes[0, 0].get_title() == "x (-0.5000 Ht)"
    assert output.is_file()
    plt.close(figure)


def test_fft_node_at_cell_centre_is_plotted_at_zero():
    field = np.zeros((30, 30))
    field[15, 15] = 1.0

    plotted, coordinates = ContourPlotter._prepare_central_plane(
        field, sz=5.0, cell_size=10.0
    )

    maximum = np.unravel_index(np.argmax(plotted), plotted.shape)
    assert coordinates[maximum[0]] == pytest.approx(0.0)
    assert coordinates[maximum[1]] == pytest.approx(0.0)
    assert coordinates[[0, -1]] == pytest.approx([-5.0, 5.0])
    assert plotted[0] == pytest.approx(plotted[-1])
    assert plotted[:, 0] == pytest.approx(plotted[:, -1])


def test_imaginary_norm_can_be_hidden(contour_data):
    system, _ = contour_data
    coordinates = np.linspace(-1.0, 1.0, 8, endpoint=False)
    z, y, x = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    result = {"orbc_rs": np.asarray([x + 1j * y + 0.1]), "orbe": [-0.2]}

    visible_figure, visible_axes = ContourPlotter.build_contourplot(
        result, system, plane="xy"
    )
    hidden_figure, hidden_axes = ContourPlotter.build_contourplot(
        result, system, plane="xy", show_imaginary_norm=False
    )

    assert "Im norm:" in visible_axes[0, 0].get_title()
    assert "Im norm:" not in hidden_axes[0, 0].get_title()
    plt.close(visible_figure)
    plt.close(hidden_figure)


def test_auto_plane_avoids_a_nodal_central_slice(contour_data):
    system, _ = contour_data
    coordinates = np.linspace(-1.0, 1.0, 8, endpoint=False)
    z, y, x = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    result = {"orbc_rs": np.asarray([z]), "orbe": [-0.2]}

    figure, axes = ContourPlotter.build_contourplot(
        result, system, plane="auto", tick_rotation=45
    )

    assert axes[0, 0].get_xlabel() == "$x$ [a.u.]"
    assert axes[0, 0].get_ylabel() == "$z$ [a.u.]"
    assert axes[0, 0].get_xticklabels()[0].get_rotation() == 45
    plt.close(figure)


@pytest.mark.parametrize("plane", ["xx", ["xy"]])
def test_contour_plotter_rejects_invalid_planes(contour_data, plane):
    system, result = contour_data
    with pytest.raises(ValueError, match="plane"):
        ContourPlotter.build_contourplot(
            result, system, plane=plane, nrows=1, ncols=2
        )
