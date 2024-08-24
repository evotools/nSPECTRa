# Changelog
All new changes are documented here.

## [v1.1.0]

### Added
- Create a "smile plot" for the derived allele frequencies (DAF)
- Demo data to test the workflow; in the future, the workflow will have CI testing to ensure basic features' stability
- Add a filtering step removing variants with derived allele frequency above a given threshold (default 0.98)

### Changed
- Docker dependencies
- DAF are computed in preprocessing, and saved as output file
- Introduced three separate options to trigger the different components: `--relate`, `--mutyper` and `--sdm` (runs `mutyper` only as default)
- Faster preprocessing of the VCF by processing by contig wherever possible (everywhere in preprocessing stage, WiP for the other components)
- Separate the filtered SDM sites based on whether they fall into a repeat masked region or not
- Increased threads provided to selected `bcftools` processes
- Ancestral genome now uses `cactus` official image, rather than on the downloaded tools
- Greatly increase performances of bed2vbed process by heavy usage of [polars](https://pola.rs/) dataframes, an improved logic and decreased I/O operations (see table):

| Version |  Mode  | Memory (GB) | Runtime (min) | Fold improv. |
|---------|--------|-------------|---------------|--------------|
| v1.0.0  |    -   |      0.7    |     18.9      |     1.0x     |
| v1.1.0  | Memory |     11.3    |      5.9      |     3.2x     |
| v1.1.0  | Greedy |     25.0    |      4.0      |     4.7x     |

### Fixed
- Fixes a start-up issue when running from the directory, rather than from `main.nf`
- Fixes `-stub` runtime issues

### Removed
- Legacy DSL1 code (`legacy/main.nf`)
- Dropped filtering process as unused in the analyses
- `BED2VBED` and `COMBINE` scripts

## [v1.0.0]
### Added
- Initial workflow release