import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pypwdft

from pypwdft import (
    DFTResult,
    GTH,
    PWDFT,
    SCFSettings,
    Structure,
)


def solver_result(norb=1, npts=5):
    coordinates = np.linspace(-1, 1, npts)
    z, y, x = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    orbitals = np.asarray([x + y + z + index for index in range(norb)], complex)
    return {
        "Etot": -1.25,
        "Ekin": 0.75,
        "Enuc": -2.0,
        "Enonloc": -0.2,
        "Erep": 0.4,
        "Exc": -0.3,
        "Eewald": 0.1,
        "edens": np.abs(orbitals[0]) ** 2 + 0.1,
        "orbe": np.linspace(-0.5, 0.1, norb),
        "orbc_rs": orbitals,
        "orbc_fft": np.fft.fftn(orbitals, axes=(1, 2, 3)),
        "converged": True,
        "iterations": 7,
        "energy_residual": 2e-7,
        "density_residual": 3e-7,
        "ttime": 1.5,
        "fft": "numpy",
        "ecut": 2.0,
        "density_ecut": 8.0,
        "npw": 33,
        "wavefunction_npts": npts,
        "density_npts": npts,
    }


def test_structure_from_name_separates_geometry_from_cutoff():
    structure = Structure.from_name("H2", cell=10)

    assert structure.symbols == ("H", "H")
    assert structure.cell == pytest.approx(10.0)
    assert structure.positions.mean(axis=0) == pytest.approx([5.0, 5.0, 5.0])
    assert "atoms=2" in repr(structure)

def test_structure_from_xyz_treats_cell_as_bohr(tmp_path):
    xyz = tmp_path / "h.xyz"
    xyz.write_text("1\nHydrogen\nH 0 0 0\n", encoding="utf-8")

    structure = Structure.from_xyz(xyz, cell=8)

    assert structure.cell == 8
    assert structure.positions[0] == pytest.approx([4, 4, 4])


def test_pwdft_builds_matching_gth_and_translates_run_options(monkeypatch):
    calculation = PWDFT(
        Structure.from_name("H2", cell=8), cutoff=2, xc="pbe"
    )
    captured = {}

    def fake_scf(**kwargs):
        captured.update(kwargs)
        return solver_result()

    monkeypatch.setattr(calculation._engine, "scf", fake_scf)
    result = calculation.run(
        convergence=2e-6,
        density_convergence=3e-6,
        max_iterations=20,
        bands=1,
        verbosity=1,
    )

    assert calculation.ionic_model == "gth-pbe"
    assert captured == {
        "tol": 2e-6,
        "density_tol": 3e-6,
        "maxiter": 20,
        "nsol": 1,
        "verbose": True,
    }
    assert isinstance(result, DFTResult)
    assert result.energy.total == -1.25
    assert result.energy.ionic == pytest.approx(-2.1)
    assert result.orbitals.energies == pytest.approx([-0.5])
    assert result.basis.plane_waves == 33
    assert result.scf.energy_residual == 2e-7
    assert result.converged
    assert result.iterations == 7
    with pytest.raises(TypeError):
        result["Etot"]
    assert not hasattr(calculation, "scf")


def test_settings_can_be_passed_to_run(monkeypatch):
    calculation = PWDFT(
        Structure.from_name("He", cell=8), cutoff=2, pseudopotential="gth"
    )
    captured = {}

    def fake_scf(**kwargs):
        captured.update(kwargs)
        return solver_result()

    monkeypatch.setattr(calculation._engine, "scf", fake_scf)
    calculation.run(SCFSettings(convergence=1e-4, max_iterations=12))

    assert captured["tol"] == 1e-4
    assert captured["maxiter"] == 12


def test_result_owns_the_common_plotting_workflow(tmp_path):
    calculation = PWDFT(Structure.from_name("He", cell=8), cutoff=2)
    result = DFTResult(
        solver_result(norb=3), system=calculation._system, noccupied=2
    )

    figure, axes = result.plot_orbitals(
        occupied=True,
        plane="xz",
        columns=2,
        save=tmp_path / "orbitals.png",
    )

    assert axes.shape == (1, 2)
    assert axes[0, 0].get_xlabel() == "$x$ [a.u.]"
    assert (tmp_path / "orbitals.png").is_file()
    plt.close(figure)


def test_advanced_and_invalid_configuration():
    structure = Structure.from_name("H2", cell=8)
    calculation = PWDFT(
        structure,
        cutoff=2,
        xc="lda",
        pseudopotential=GTH(charges={"H": 1}),
        fft_backend="numpy",
    )
    assert calculation.ionic_model == "gth-pade"

    with pytest.raises(ValueError, match="does not match"):
        PWDFT(structure, cutoff=2, xc="pbe", pseudopotential="gth-pade")
    with pytest.raises(ValueError, match="requires device='cuda'"):
        PWDFT(structure, cutoff=2, fft_backend="cupy")


def test_legacy_api_is_not_publicly_exported():
    for name in (
        "PyPWDFT",
        "PeriodicSystem",
        "SystemBuilder",
        "GTHPseudopotential",
        "ContourPlotter",
        "plot_orbital_contours",
    ):
        assert not hasattr(pypwdft, name)
