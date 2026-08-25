# Contributing

Changes should be made through a branch and pull request. Keep biological method changes separate from mechanical refactoring and document any altered coordinate, strand, counting, filtering or statistical contract.

Before requesting review:

```bash
bash tests/run_tests.sh
python3 -m build
```

New protocol profiles require synthetic coordinate/orientation truth tests. New reference profiles require assembly and contig validation. Changes to APA-B must preserve independence from APA-A and repeat the no-UMI/no-dedup REV pilot.

Never commit sequencing data, alignments, reference genomes, annotations, PAS atlases, model weights, credentials, personal paths or project results.
