# Related work: tumor-informed ctDNA MRD landscape

nf-mrd-umi is a two-pipeline, tumor-informed MRD design. Pipeline A discovers a
patient-specific somatic mutation set from tumor+normal WES and designs a locked
interrogation panel; Pipeline B takes cfDNA through duplex-UMI consensus, then
interrogates (does NOT re-call) that locked panel and emits a panel-integrated
MRD verdict using an empirical Monte-Carlo null plus an odds-space,
allele-specific-enrichment de-biasing step. nf-mrd-umi deliberately does NOT
re-implement standard genomics: it delegates alignment, somatic calling, and
UMI consensus to the same biocontainers that nf-core/sarek and nf-core/fastquorum
use, and only adds the MRD-specific layers (locked-panel interrogation, empirical
null, enrichment de-biasing) in bin/. This document maps each component to the
closest verified open-source building block, and contrasts the design with the
commercial CLIA MRD services for context. All repos below were resolved by
fetching their github.com pages (or, where noted, by web-search listing).

## Commercial tumor-informed MRD assays

These are CLIA laboratory SERVICES, not downloadable software. Listed for context
on assay strategy only.

| Assay | Vendor | Approach | Allele-specific enrichment? | Software available? |
|-------|--------|----------|-----------------------------|---------------------|
| Signatera | Natera | Tumor-informed bespoke mPCR (~16 variants) + deep NGS | No (mPCR amplicon deep-seq) | No, CLIA service |
| RaDaR (ex-Inivata) | NeoGenomics | Tumor-informed personalized panel (up to ~48 variants), deep NGS | No (anchored-multiplex deep-seq) | No, CLIA service |
| NeXT Personal | Personalis | Tumor-informed ultra-broad WGS-derived panel (~1.8k-25k variants), deep NGS | No (deep-seq of large panel) | No, CLIA service |
| FoundationOne Tracker | Foundation Medicine | Tumor-informed; tissue CDx selects variants, monitored by deep NGS (mPCR) | No (deep-seq) | No, CLIA service |
| Haystack MRD | Haystack Oncology (Quest) | Tumor-informed bespoke panel, duplex deep-seq | No (duplex deep-seq) | No, CLIA service |
| C2i | C2i Genomics | Tumor-informed whole-genome aggregation (MRDetect-style) | No (WGS signal integration) | No, CLIA service |
| Reveal | Guardant | Tumor-naive, blood-only; genomic + methylation/epigenomic | No (methylation/deep-seq) | No, CLIA service |
| xM | Tempus | Tumor-naive, blood-only; methylation + variant MRD | No (methylation/deep-seq) | No, CLIA service |

Note: essentially NONE of these use allele-specific enrichment chemistry; they are
deep-seq of a personalized panel (mPCR/hybrid-capture/duplex) or genome-wide WGS
aggregation, or (Guardant/Tempus) tumor-naive methylation. Allele-specific
enrichment of the mutant allele is closer to NaME-PrO (Song et al., NAR 2016) and
COLD-PCR (Li et al., Nat Med 2008) chemistries, which are uncommon in the MRD
assay market. The nf-mrd-umi enrichment de-biasing layer is built to model exactly
that kind of allele-specific over-representation, which is what makes it distinct.

## Open-source building blocks (by pipeline component)

Mapping nf-mrd-umi parts to the closest verified open repo.

### Somatic discovery (Pipeline A)
| Component | Repo | URL | Lang / License | What it does |
|-----------|------|-----|----------------|--------------|
| End-to-end somatic pipeline | nf-core/sarek | https://github.com/nf-core/sarek | Nextflow / MIT | Germline+somatic variant calling from WGS/WES/targeted (pre-proc, calling, annotation) |
| Somatic SNV/indel caller | broadinstitute/gatk (Mutect2) | https://github.com/broadinstitute/gatk | Java / Apache-2.0 | GATK4; Mutect2 is the somatic caller used for tumor/normal discovery |

### Duplex / UMI consensus (Pipeline B)
| Component | Repo | URL | Lang / License | What it does |
|-----------|------|-----|----------------|--------------|
| UMI consensus pipeline | nf-core/fastquorum | https://github.com/nf-core/fastquorum | Nextflow / MIT | Produces consensus reads from UMIs (wraps fgbio best-practice duplex workflow) |
| UMI consensus toolkit | fulcrumgenomics/fgbio | https://github.com/fulcrumgenomics/fgbio | Scala / MIT | GroupReadsByUmi, Call(Duplex)ConsensusReads -- the core duplex consensus engine |
| UMI handling | CGATOxford/UMI-tools | https://github.com/CGATOxford/UMI-tools | Python / MIT | UMI extraction/dedup/grouping (directional method) |
| UMI error-correct + call | stahlberggroup/umierrorcorrect | https://github.com/stahlberggroup/umierrorcorrect | Python / open (see repo) | fastq-to-consensus UMI pipeline with error correction and calling |
| UMI assemble+align+call | mikessh/mageri | https://github.com/mikessh/mageri | Java / see repo | MAGERI: molecular-barcode consensus assembly, alignment and variant calling |

### Error model / empirical null / low-VAF callers
| Component | Repo | URL | Lang / License | What it does |
|-----------|------|-----|----------------|--------------|
| Read-level error model | JakobSkouPedersenLab/dreams | https://github.com/JakobSkouPedersenLab/dreams | R / GPL-3.0 | DREAMS: deep read-level error model for ctDNA, closest published analogue to an empirical null |
| Beta-binomial subclonal caller | gerstung-lab/deepSNV | https://github.com/gerstung-lab/deepSNV | R/C / GPL (see repo) | deepSNV/shearwater: beta-binomial LRT vs a control panel for low-VAF SNVs |
| UMI low-VAF caller | (smCounter2) qiaseq/smcounter-v2-paper | https://github.com/qiaseq/smcounter-v2-paper | Python / see repo | Paper code for smCounter2 (beta-binomial UMI caller). NOTE: smCounter2 itself (qiaseq/qiaseq-dna) is no longer open-source |
| Original smCounter | xuchang116/smCounter | https://github.com/xuchang116/smCounter | Python / see repo | smCounter v1 UMI-aware low-VAF caller (predecessor) |
| Poisson ctDNA SNV caller | sfu-compbio/sinvict | https://github.com/sfu-compbio/sinvict | C++ / see repo | SiNVICT: ultra-sensitive SNV/indel detection in cfDNA at low VAF |
| Amplicon/UMI caller | AstraZeneca-NGS/VarDictJava | https://github.com/AstraZeneca-NGS/VarDictJava | Java / MIT | VarDict low-frequency caller (repo archived May 2026) |

### Tumor-fraction / WGS-MRD
| Component | Repo | URL | Lang / License | What it does |
|-----------|------|-----|----------------|--------------|
| Tumor fraction from ULP-WGS | GavinHaLab/ichorCNA | https://github.com/GavinHaLab/ichorCNA | R / GPL-3.0 | Estimates tumor fraction from ultra-low-pass WGS (fork of broadinstitute/ichorCNA) |
| WGS tumor-informed MRD | oicr-gsi/mrdetect | https://github.com/oicr-gsi/mrdetect | workflow / see repo | OICR workflow wrapping MRDetect (Landau-lab genome-wide cfDNA integration). Verified via web-search listing; page fetch intermittently 404s |

### CNV / clonality (Pipeline A)
| Component | Repo | URL | Lang / License | What it does |
|-----------|------|-----|----------------|--------------|
| Allele-specific CN + purity | mskcc/facets | https://github.com/mskcc/facets | R/C / see repo | FACETS: fraction and copy-number estimate from tumor/normal seq |
| Clonal structure | Roth-Lab/pyclone-vi | https://github.com/Roth-Lab/pyclone-vi | Python / GPL-3.0 | PyClone-VI: fast clonal-population inference from SNVs (clusters variants for panel selection) |

Not separately packaged as a public repo: CAPP-Seq / iDES (Stanford method,
commercialized as Roche AVENIO -- no canonical open repo) and MRD-EDGE (Landau lab
successor to MRDetect -- published in Nature Medicine 2024 but code not publicly
released as of this writing). Treat both as "unverified / no open repo".

## The gap nf-mrd-umi fills

No open repo packages the full nf-mrd-umi recipe: tumor-informed interrogation at a
patient-LOCKED panel + panel-integrated Monte-Carlo empirical null +
odds-space allele-specific-enrichment de-biasing, as a turnkey pipeline. The
individual pieces exist (sarek/Mutect2 for discovery; fastquorum/fgbio for duplex
consensus; DREAMS/deepSNV/SiNVICT for low-VAF error modeling; ichorCNA/MRDetect for
WGS tumor fraction; FACETS/PyClone-VI for CN and clonality), but they are
single-purpose tools, not an integrated tumor-informed MRD caller. The realistic
open recipe is: sarek (A) + fastquorum (B) + a custom MRD caller -- which is exactly
what nf-mrd-umi's bin/ provides. Crucially, the existing low-VAF callers all assume
a fixed, unbiased error background; none model allele-specific over-representation
of the mutant allele. The odds-space allele-specific-enrichment de-biasing layer is
the genuinely distinctive piece, and it has no open-source equivalent in the MRD
tooling space (its conceptual kin is NaME-PrO / COLD-PCR enrichment chemistry, not
any MRD software).

## Links

- nf-core/sarek -- https://github.com/nf-core/sarek
- broadinstitute/gatk (Mutect2) -- https://github.com/broadinstitute/gatk
- nf-core/fastquorum -- https://github.com/nf-core/fastquorum
- fulcrumgenomics/fgbio -- https://github.com/fulcrumgenomics/fgbio
- CGATOxford/UMI-tools -- https://github.com/CGATOxford/UMI-tools
- stahlberggroup/umierrorcorrect -- https://github.com/stahlberggroup/umierrorcorrect
- mikessh/mageri -- https://github.com/mikessh/mageri
- JakobSkouPedersenLab/dreams -- https://github.com/JakobSkouPedersenLab/dreams
- gerstung-lab/deepSNV -- https://github.com/gerstung-lab/deepSNV
- qiaseq/smcounter-v2-paper -- https://github.com/qiaseq/smcounter-v2-paper
- xuchang116/smCounter -- https://github.com/xuchang116/smCounter
- sfu-compbio/sinvict -- https://github.com/sfu-compbio/sinvict
- AstraZeneca-NGS/VarDictJava -- https://github.com/AstraZeneca-NGS/VarDictJava
- GavinHaLab/ichorCNA -- https://github.com/GavinHaLab/ichorCNA
- oicr-gsi/mrdetect -- https://github.com/oicr-gsi/mrdetect
- mskcc/facets -- https://github.com/mskcc/facets
- Roth-Lab/pyclone-vi -- https://github.com/Roth-Lab/pyclone-vi
