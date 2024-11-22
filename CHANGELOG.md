# Changelog
All new changes are documented here.

## [v1.1.1]
### Changed
- The workflow now filters out all sites where the ancestral allele does not match either REF or ALT in the VCF
- Heavy code maintainance, with much better code linting

### Fixed
- Few cases where the workflow would mix up the DAF

## [v1.1.0]

### Added
- Create a "smile plot" for the derived allele frequencies (DAF)
- Demo data to test the workflow; in the future, the workflow will have CI testing to ensure basic features' stability
- Add a filtering step removing variants with derived allele frequency above a given threshold (default 0.98)

### Changed
- Docker dependencies
- DAF are computed in preprocessing, and saved as output file
- Introduced three separate options to trigger the different components: `--relate`, `--mutyper` and `--sdm` (runs `mutyper` only as default)
- `--relate_path` is now used to provide the path to the relate installation directory instead of `--relate`
- `--ancestral` is now `--ancestral_fna`, and `--ref_fasta` is now `--fasta_fna`
- Faster preprocessing of the VCF by processing by contig wherever possible
- Faster VCF I/O thanks to dropping most `INFO` fields when extracting biallelic sites
- Separate the filtered SDM sites based on whether they fall into a repeat masked region or not
- Increased threads provided to selected `bcftools` processes
- Ancestral genome now uses `cactus` official image, rather than on the downloaded tools
- The workflow now uses chunking whereever possible to speed up processing, defined by `chunk_size`
    - The chunk size can be slightly higher when consecutive sites are found, effectively splitting only when breaks in variants are identified
    - Both Mutyper and SDM subworkflows takes advantage of the approach, reducing redundancy of the analyses and allowing lower I/O with faster analysis of data
- Collection of mutation type is now performed by sequence, rather on the full dataset  
- Greatly increase performances of bed2vbed process by heavy usage of [polars](https://pola.rs/) dataframes, an improved logic and decreased I/O operations (see table):

| Version |  Mode  | Memory (GB) | Runtime (min) | Fold improv. |
|---------|--------|-------------|---------------|--------------|
| v1.0.0  |    -   |      0.7    |     18.9      |     1.0x     |
| v1.1.0  | Memory |     11.3    |      5.9      |     3.2x     |
| v1.1.0  | Greedy |     25.0    |      4.0      |     4.7x     |

### Fixed
- Fixes a start-up issue when running from the directory, rather than from `main.nf`
- Fixes `-stub` runtime issues
- Mutyper plots ignoring groupings

### Removed
- Legacy DSL1 code (`legacy/main.nf`)
- Dropped filtering process as unused in the analyses
- `BED2VBED` and `COMBINE` scripts

## [v1.0.0]
### Added
- Initial workflow release