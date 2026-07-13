# Bundled GTH pseudopotentials

The `pade` and `pbe` directories contain the LDA/PADE and PBE
Goedecker-Teter-Hutter pseudopotentials distributed with eminus 3.2.2. The
parameter files originate from the CP2K data repository. See `NOTICE` for
provenance and licensing.

Selecting `pseudopotential="gth"` in `PWDFT` automatically chooses the family
matching the exchange-correlation functional. Use `GTH(family="pade")` or
`GTH(family="pbe")` for an explicit selection. `GTH(path=...)` reads another
CP2K-style parameter directory.
