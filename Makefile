SHELL := /bin/bash
PACKAGE=depot
VERSION=$(shell grep Version pkg/DESCRIPTION |awk '{print $$2}')
TARBALL=${PACKAGE}_${VERSION}.tar.gz
PKGDIR=.

all:
	Rscript -e 'devtools::document("pkg")'
	R CMD build pkg

bump:
	Rscript pkg/inst/scripts/release.R

release:
	make all
	cp ${TARBALL} src/contrib
	Rscript -e 'tools::write_PACKAGES("src/contrib/")'
