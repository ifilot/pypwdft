import matplotlib.pyplot as plt

from pypwdft import PWDFT, Structure


def main():
    structure = Structure.from_name("benzene", cell=15)
    calculation = PWDFT(
        structure,
        cutoff=40,              # wavefunction cutoff in Hartree
        xc="pbe",               # exchange-correlation functional
        pseudopotential="gth",  # pseudopotential family
        device="cuda",          # device to use for calculations (cpu or cuda)
    )
    result = calculation.run(
        convergence=1e-4,           # convergence threshold for self-consistent field iterations
        density_convergence=1e-4,   # convergence threshold for density matrix updates
        verbosity=1,                # verbosity level for output (0: silent, 1: normal, 2: verbose)
        bands=18                    # number of bands to compute (including unoccupied orbitals)
    )

    print("Total energy:", result.energy.total)
    result.plot_orbitals(
        plane="auto",               # avoid plotting an orbital's nodal plane
        save="examples/example_benzene_contours.png",
        show_imaginary_norm=False,
        columns=6,
        ngrid=3,
        tick_rotation=45,
        occupied=False,             # also plot unoccupied orbitals
    )
    plt.show()


if __name__ == "__main__":
    main()
