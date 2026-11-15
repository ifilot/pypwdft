from pathlib import Path

import matplotlib.pyplot as plt

from pypwdft import PWDFT, Structure


def main():
    structure = Structure.from_name("CO", cell=10)
    calculation = PWDFT(
        structure,
        cutoff=40,              # wavefunction cutoff in Hartree
        xc="pbe",               # exchange-correlation functional
        pseudopotential="gth",  # pseudopotential family
        device="cpu",           # "cuda" runs on a GPU (requires CuPy)
    )
    result = calculation.run(
        convergence=1e-4,           # SCF total-energy convergence threshold (Ha)
        density_convergence=1e-4,   # SCF RMS density-residual threshold
        verbosity=1,                # print SCF progress (0: silent)
        bands=6
    )

    print("Total energy:", result.energy.total)
    result.plot_orbitals(
        plane="xz",
        save=Path(__file__).with_name("example_co_contours.png"),
        show_imaginary_norm=False,
        occupied=False,
    )
    plt.show()


if __name__ == "__main__":
    main()
