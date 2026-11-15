"""Benchmark the SCF solver across eigensolvers, devices and FFT threads.

Every configuration is run ``--repeats`` times and the best wall time is
reported, together with two hardware-independent figures: the number of SCF
iterations and the number of Hamiltonian applications.

Examples
--------
CPU, both eigensolvers, one and eight pyFFTW threads::

    python benchmarks/scf_benchmark.py --threads 1 8

GPU (requires CuPy and a CUDA device)::

    python benchmarks/scf_benchmark.py --device cuda

A heavier case::

    python benchmarks/scf_benchmark.py --molecule benzene --cell 16 \\
        --cutoff 30 --bands 18 --device cuda
"""

import argparse
import time

from pypwdft import PWDFT, Structure


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
    )
    parser.add_argument("--molecule", default="CO",
                        help="bundled structure name (default: CO)")
    parser.add_argument("--cell", type=float, default=10.0,
                        help="cubic cell edge in bohr (default: 10)")
    parser.add_argument("--cutoff", type=float, default=40.0,
                        help="wavefunction cutoff in Hartree (default: 40)")
    parser.add_argument("--xc", default="pbe", choices=("lda", "pbe"))
    parser.add_argument("--bands", type=int, default=None,
                        help="number of orbitals (default: occupied only)")
    parser.add_argument("--convergence", type=float, default=1e-6,
                        help="energy and density threshold (default: 1e-6)")
    parser.add_argument("--device", default="cpu", choices=("cpu", "cuda"))
    parser.add_argument("--fft-backend", default=None,
                        choices=("pyfftw", "numpy", "scipy", "cupy"),
                        help="override the device's default FFT backend")
    parser.add_argument("--eigensolvers", nargs="+",
                        default=("davidson", "arpack"),
                        choices=("davidson", "arpack"))
    parser.add_argument("--threads", nargs="+", type=int, default=(1,),
                        help="pyFFTW thread counts to compare (CPU only)")
    parser.add_argument("--repeats", type=int, default=3,
                        help="runs per configuration; the best is reported")
    return parser.parse_args()


def main():
    args = parse_arguments()
    structure = Structure.from_name(args.molecule, cell=args.cell)
    uses_fftw = args.device == "cpu" and args.fft_backend in (None, "pyfftw")
    thread_counts = args.threads if uses_fftw else (1,)

    rows = []
    for eigensolver in args.eigensolvers:
        for threads in thread_counts:
            calculation = PWDFT(
                structure,
                cutoff=args.cutoff,
                xc=args.xc,
                device=args.device,
                fft_backend=args.fft_backend,
                threads=threads,
            )
            times = []
            for _ in range(args.repeats):
                start = time.perf_counter()
                result = calculation.run(
                    convergence=args.convergence,
                    bands=args.bands,
                    eigensolver=eigensolver,
                )
                times.append(time.perf_counter() - start)
            rows.append((eigensolver, threads, min(times), result))
            print(
                f"finished {eigensolver}, {threads} thread(s): "
                + ", ".join(f"{value:.2f} s" for value in times),
                flush=True,
            )

    basis = rows[0][3].basis
    print()
    print(
        f"{args.molecule}, {args.xc.upper()}, cell {args.cell:g} bohr, "
        f"cutoff {args.cutoff:g} Ha, convergence {args.convergence:g}"
    )
    print(
        f"backend {calculation.fft_backend}, grid {basis.density_grid}^3, "
        f"{basis.plane_waves} plane waves, "
        f"{len(rows[0][3].orbitals.energies)} orbitals"
    )
    print()
    print("eigensolver  threads  best time (s)  SCF iter  H applications"
          "     total energy (Ha)")
    for eigensolver, threads, best, result in rows:
        print(
            f"{eigensolver:<11s}  {threads:7d}  {best:13.2f}  "
            f"{result.iterations:8d}  "
            f"{result.scf.hamiltonian_applications:14d}  "
            f"{result.energy.total:20.10f}"
        )


if __name__ == "__main__":
    main()
