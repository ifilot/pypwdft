import matplotlib.pyplot as plt
from pypwdft import PWDFT, Structure

# Purpose: perform a self-consistent field calculation of H2 molecule and
# visualize its valence orbitals. This script uses the PWDFT class to perform
# the calculation and generate contour plots of the molecular orbitals.


def main():
    structure = Structure.from_name("H2", cell=10)
    calculation = PWDFT(
        structure,
        cutoff=10,              # wavefunction cutoff in Hartree
        xc="pbe",               # exchange-correlation functional
        pseudopotential="gth",  # pseudopotential family
        device="cpu",           # device to use for calculations (cpu or cuda)
    )
    result = calculation.run(
        convergence=1e-6,
        density_convergence=1e-6,
        verbosity=1,
        bands=2,
    )

    print("Total energy:", result.energy.total)
    result.plot_orbitals(
        plane="xz",
        save="examples/example_h2_contours.png",
        occupied=False,
    )
    plt.show()


if __name__ == "__main__":
    main()
