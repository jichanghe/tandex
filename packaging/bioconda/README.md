# Bioconda recipe

Not yet submitted. To publish tandex on Bioconda:

1. Tag and publish the release on GitHub (for example `v0.7.0`).
2. Put the release tarball's sha256 into `meta.yaml` (the command is in the file).
3. Fork https://github.com/bioconda/bioconda-recipes, copy `meta.yaml` and `build.sh` to
   `recipes/tandex/`, and open a pull request. Bioconda's own checks build and test the package.

The recipe has not been built with conda yet; the build command in `build.sh` is the same one
the Makefile uses.
