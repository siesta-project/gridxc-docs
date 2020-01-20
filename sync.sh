#!/bin/sh
#
# Use this file to sync the contents of FORD's target directory
# and the sources for the GitHub Pages documentation
#
# Note the exclusions, to avoid removing the Git control file and this very same file.
#
FORD_TARGET=/tmp/gridxc-docs/   # Note trailing slash
rsync -av --delete --exclude=README.md \
          --exclude=.git \
          --exclude=sync.sh ${FORD_TARGET} .
