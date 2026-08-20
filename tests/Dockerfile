# Test image for the tool-level characterization tests.
#
#   docker build --platform linux/amd64 -t wurare-svdb-tests tests/
#   docker run --rm --platform linux/amd64 -v "$PWD:/w" -w /w wurare-svdb-tests \
#       python3 tests/test_svdb_assumptions.py
#
# Why an image rather than the conda env in environment.yml: bioconda ships svdb for
# linux-64 only, so `conda env create` cannot produce a working svdb on an arm64 Mac. The
# pipeline runs from containers on the cluster anyway, which is what environment.yml's own
# header already says about svtk, duphold, truvari and AnnotSV -- svdb joins that list on
# this architecture. --platform linux/amd64 is required and emulated; the tests are small.
#
# Why not a pinned biocontainer: the svdb one carries svdb without bcftools and the bcftools
# one the reverse, and test_svdb_assumptions.py needs both on one PATH to read a merged VCF
# back. Versions here track environment.yml, not the pipeline's per-process containers --
# this image tests the *tool*, and the pipeline pins its own.
FROM mambaorg/micromamba:1.5.8

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        svdb=2.12.0 \
        bcftools=1.24 \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
