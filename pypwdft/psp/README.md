# Bundled GTH pseudopotentials

The `pade` and `pbe` directories contain the LDA/PADE and PBE
Goedecker-Teter-Hutter pseudopotentials distributed with eminus 3.2.2. The
parameter files originate from the CP2K data repository. See `NOTICE` for
provenance and licensing.

Select them with `GTHPseudopotential(family="pade")` or
`GTHPseudopotential(family="pbe")`. The former remains the default.
`GTHPseudopotential(path=...)` can read another CP2K-style GTH parameter
directory explicitly.
