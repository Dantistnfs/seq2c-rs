
Rewrite of [seq2c's](https://github.com/AstraZeneca-NGS/Seq2C) seq2cov.pl script in rust.


Efficently uses parallel cpus and greatly improves performance for big .bed panels. (35 seconds vs 90 minutes)


## Usage

Mimics seq2cov 1-1 in terms of input and output for better compatibility.
```bash 
seq2c-rs -b path_to_bam/sample.bam -N sample_name -p panel.bed --threads 16 > output.tsv
```
## Benchmark

Bam file ~15Gb

| Bed file (number of rows) | Perl (requires index) | Rust        | Rust v1.0.0 threads=16 | Rust v1.0.1 threads=16 |
| ------------------------- | --------------------- | ----------- | ---------------------- | ---------------------- |
| 10                        | 0.409s                | 270s (4.5m) | ~74 seconds            | 34.95s                 |
| 100                       | 3.058s                | 270s (4.5m) | ~74 seconds            | 33.93s                 |
| 1000                      | 23.181s               | 270s (4.5m) | ~74 seconds            | 31.65s                 |
| 10000                     | ~240s (4 min)         | 270s (4.5m) | ~74 seconds            | 33.0s                  |
| 100000                    | ~40m                  | 270s (4.5m) | ~74 seconds            | 34.486s                |
| 261643                    | ~90m                  | 270s (4.5m) | 74 seconds             | 35.9s                  |

