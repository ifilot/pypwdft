"""Study plane-wave cutoff convergence for CH4 with GTH-PBE.

Each completed cutoff is written immediately to a CSV file.  Re-running this
script therefore resumes the sweep instead of repeating finished calculations.
"""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pypwdft import PWDFT, Structure


HARTREE_TO_EV = 27.211386245988
CACHE_FIELDS = (
    "ecut_ha",
    "ecut_ev",
    "energy_ha",
    "npw",
    "wavefunction_npts",
    "density_npts",
    "density_ecut_ha",
    "iterations",
    "energy_residual_ha",
    "density_residual",
    "time_s",
    "backend",
)


def read_cache(path):
    """Read cached results, indexed by their integer cutoff."""
    if not path.exists():
        return {}

    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        missing = set(CACHE_FIELDS) - set(rows.fieldnames or ())
        if missing:
            raise ValueError(
                f"Cache file {path} is missing columns: {sorted(missing)}"
            )
        return {int(float(row["ecut_ha"])): row for row in rows}


def write_cache(path, results):
    """Atomically replace the cache with the currently completed results."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CACHE_FIELDS)
        writer.writeheader()
        for ecut in sorted(results):
            writer.writerow(results[ecut])
    temporary.replace(path)


def calculate_point(molecule, ecut, backend):
    """Run one GTH-PBE/PBE calculation and return its scalar results."""
    structure = Structure.from_name(molecule, cell=10)
    calculation = PWDFT(
        structure,
        cutoff=ecut,
        xc="pbe",
        pseudopotential="gth",
        device="cuda" if backend == "cupy" else "cpu",
        fft_backend=backend,
    )
    result = calculation.run(
        convergence=1e-5,
        density_convergence=1e-5,
        max_iterations=100,
        verbosity=1,
    )
    return {
        "ecut_ha": ecut,
        "ecut_ev": f"{ecut * HARTREE_TO_EV:.8f}",
        "energy_ha": f"{result.energy.total:.12f}",
        "npw": result.basis.plane_waves,
        "wavefunction_npts": result.basis.wavefunction_grid,
        "density_npts": result.basis.density_grid,
        "density_ecut_ha": f"{result.basis.density_cutoff:.8f}",
        "iterations": result.scf.iterations,
        "energy_residual_ha": f"{result.scf.energy_residual:.12e}",
        "density_residual": f"{result.scf.density_residual:.12e}",
        "time_s": f"{result.scf.elapsed_time:.6f}",
        "backend": backend,
    }


def print_results(molecule, results):
    """Print a compact convergence table."""
    print(f"\n{molecule} cutoff-convergence summary")
    print("-" * 102)
    print(
        " Ecut (Ha)  Ecut (eV)       Npw   WF grid  Density grid"
        "       Energy (Ha)   dE from previous (mHa)"
    )
    print("-" * 102)
    previous = None
    for ecut in sorted(results):
        row = results[ecut]
        energy = float(row["energy_ha"])
        difference = "--" if previous is None else f"{(energy - previous) * 1e3:10.4f}"
        print(
            f"{ecut:10d}  {float(row['ecut_ev']):9.2f}  "
            f"{int(row['npw']):9d}  {int(row['wavefunction_npts']):5d}^3  "
            f"{int(row['density_npts']):10d}^3  {energy:16.10f}  "
            f"{difference:>22s}"
        )
        previous = energy
    print("-" * 102)


def plot_results(molecule, results, path):
    """Plot the energy convergence and plane-wave basis growth."""
    cutoffs = np.array(sorted(results), dtype=float)
    energies = np.array(
        [float(results[int(ecut)]["energy_ha"]) for ecut in cutoffs]
    )
    plane_waves = np.array(
        [int(results[int(ecut)]["npw"]) for ecut in cutoffs]
    )
    error_ev = np.abs(energies - energies[-1]) * HARTREE_TO_EV

    figure, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    axes[0].plot(cutoffs, energies, "o-", color="tab:blue")
    axes[0].set_xlabel("Wavefunction cutoff (Ha)")
    axes[0].set_ylabel("Total energy (Ha)")
    molecule_label = molecule.replace("CH4", "CH$_4$")
    axes[0].set_title(f"{molecule_label} total energy")

    if len(cutoffs) > 1:
        axes[1].semilogy(
            cutoffs[:-1], error_ev[:-1], "o-", color="tab:orange"
        )
    axes[1].axvline(
        cutoffs[-1], color="0.5", linestyle="--", linewidth=1
    )
    axes[1].set_xlabel("Wavefunction cutoff (Ha)")
    axes[1].set_ylabel(
        f"|E - E({cutoffs[-1]:g} Ha)| (eV)"
    )
    axes[1].set_title("Energy convergence")

    axes[2].plot(cutoffs, plane_waves, "o-", color="tab:green")
    axes[2].set_xlabel("Wavefunction cutoff (Ha)")
    axes[2].set_ylabel("Number of plane waves")
    axes[2].set_title("Basis size")
    axes[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    for axis in axes:
        axis.grid(alpha=0.25)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    return figure


def parse_args(molecule):
    example_directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=(
            f"Study plane-wave cutoff convergence for {molecule} with "
            "GTH-PBE. Completed points are cached in CSV format."
        )
    )
    parser.add_argument(
        "--backend",
        choices=("cupy", "pyfftw", "numpy", "scipy"),
        default="pyfftw",
        help="FFT/linear-algebra backend (default: pyfftw).",
    )
    parser.add_argument(
        "--max-ecut",
        type=int,
        default=120,
        help="Largest cutoff in Ha; cutoffs start at 20 Ha (default: 120).",
    )
    parser.add_argument(
        "--cache",
        type=Path,
        default=(
            example_directory
            / "data"
            / f"{molecule.lower()}_ecut_convergence.csv"
        ),
        help="CSV cache path.",
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=(
            example_directory
            / "data"
            / f"{molecule.lower()}_ecut_convergence.png"
        ),
        help="Output figure path.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Save the figure without opening a plot window.",
    )
    args = parser.parse_args()
    if args.max_ecut < 20 or args.max_ecut % 20:
        parser.error("--max-ecut must be a positive multiple of 20 Ha")
    return args


def main(molecule="CH4"):
    args = parse_args(molecule)
    cutoffs = range(20, args.max_ecut + 1, 20)
    results = read_cache(args.cache)

    for ecut in cutoffs:
        if ecut in results:
            print(f"Reusing cached ecut={ecut} Ha from {args.cache}")
            continue
        print(f"\nCalculating {molecule} with ecut={ecut} Ha")
        results[ecut] = calculate_point(molecule, ecut, args.backend)
        write_cache(args.cache, results)
        print(f"Stored completed ecut={ecut} Ha point in {args.cache}")

    selected_results = {ecut: results[ecut] for ecut in cutoffs}
    print_results(molecule, selected_results)
    plot_results(molecule, selected_results, args.figure)
    print(f"Figure written to {args.figure}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
