# Containers

The container strategy mirrors [SVcaller](https://github.com/alvin8-git/SVcaller)
and exists to avoid two failure modes seen with monolithic images: **solver
lockups** and **size ballooning**.

## The principle

**One pinned public biocontainer per process. Never a single fat conda image.**

A monolithic image that conda-installs every tool re-solves a huge dependency
graph at build time (the mamba/conda solver lockup) and stacks every tool's
dependencies into one layer (the balloon). A per-tool biocontainer ships a
finished layer where that solve already happened upstream, so Nextflow just
pulls it. Conda stays available only as a `local` fallback profile, never the
default.

## How it's wired

- **`conf/docker.config`** assigns one container per process via `withName:`
  selectors, pinned by full tag including the conda build hash (verified live on
  quay.io/biocontainers). Modules carry only a resource `label`.
- **`conf/base.config`** sets the default container (the custom engine image) and
  the `process_low/medium/high` resource labels.
- **`conf/local.config`** is the conda fallback (`/data/alvin/envs/mrd`), Docker
  off.
- **`nextflow.config`** wires the `docker` / `local` / `test` profiles.

### Pinned tool images

| Process | Image |
|---|---|
| `BWAMEM2_MEM` | `bwa-mem2:2.2.1--he70b90d_8` |
| `SAMTOOLS_*` | `samtools:1.23.1--ha83d96e_0` |
| `FGBIO_*` | `fgbio:4.1.0--hdfd78af_0` |
| `GATK4_*` (Mutect2) | `gatk4:4.6.2.0--py310hdfd78af_1` |
| `PICARD_*` | `picard:3.4.0--hdfd78af_0` |
| `VEP` | `ensembl-vep:116.0--pl5321h2a3209d_0` |
| `CNVKIT` | `cnvkit:0.9.13--pyhdfd78af_0` |
| `PURPLE` | `hmftools-purple:4.4--hdfd78af_0` |
| `FACETS` | `cnv_facets:0.16.1--py312r43h8537716_1` |
| `PYCLONEVI` | `pyclone-vi:0.2.0--pyhdfd78af_0` |
| engine (`INTERROGATE`, `MRD_INTEGRATE`, `PANEL_SELECT`, `BUILD_BACKGROUND`, `BUILD_PANEL`, `VALIDATE`) | `mrd-umi/utils:1.0` (custom) |

All biocontainers prefixed `quay.io/biocontainers/`.

## The one custom image (`Dockerfile.utils`)

Built only because the MRD engine (custom `bin/` + numpy/scipy/pysam) is not a
biocontainer. Anti-balloon techniques:

- `FROM python:3.11-slim` — slim base, never a full OS or conda env.
- single `RUN apt-get install --no-install-recommends ... && rm -rf /var/lib/apt/lists/*`.
- `pip install --no-cache-dir` with every dependency pinned `==`.
- `COPY bin/` last; `.dockerignore` keeps sequencing data, references, run dirs,
  and git history out of the build context.

```bash
docker build -f Dockerfile.utils -t mrd-umi/utils:1.0 .
```

## Rule of thumb

If a step needs a tool that isn't in its container (e.g. bwa-mem2 has no
samtools), add a **separate process** with its own single-tool container. Do not
rebuild a combined image — that reintroduces the lockup and the balloon.
