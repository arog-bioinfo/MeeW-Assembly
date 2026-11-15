# Synthetic smoke-test reads

All six FASTQ files in this directory are tiny, synthetic inputs created only to exercise MeeW-Assembly workflow routing and tool invocation:

| Mode | Files | Role |
| --- | --- | --- |
| Short metaSPAdes | `short_sample_R1.fastq`, `short_sample_R2.fastq` | Paired short reads |
| Long-Nanopore Flye | `long_sample.fastq` | Nanopore-like long reads |
| Hybrid-PacBio metaSPAdes | `hybrid_sample_R1.fastq`, `hybrid_sample_R2.fastq`, `hybrid_sample_long.fastq` | Paired short reads plus PacBio-like long reads |

The reads are non-biological and must not be used for scientific interpretation, benchmarking, or tool-accuracy evaluation. They are committed so smoke tests are deterministic and network-independent: CI never downloads fixture data from public or biological datasets. For the same reason, each test configuration uses its generated synthetic assembly as a smoke-only QUAST reference; those self-reference metrics are not scientific quality measurements.
