"""Profile one SCF run: python profile_case.py <device> <molecule> <cell> <cutoff>"""
import cProfile, pstats, sys
from pypwdft import PWDFT, Structure
dev, mol, cell, cutoff = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
c = PWDFT(Structure.from_name(mol, cell=cell), cutoff=cutoff, device=dev,
          threads=8 if dev == "cpu" else None)
if dev == "cuda":
    c.run(convergence=1e-6)     # warm-up
p = cProfile.Profile(); p.enable(); c.run(convergence=1e-6); p.disable()
pstats.Stats(p).sort_stats("tottime").print_stats(14)
