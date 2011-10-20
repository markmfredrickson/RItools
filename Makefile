R: RITOOLS_TIMESTAMP
	R_PROFILE=Rprofile R -q --no-save 

.local:
	mkdir .local

RITOOLS_TIMESTAMP: .local R/* tests/* 
	R --vanilla CMD Install --library=.local .
	date > RITOOLS_TIMESTAMP

autotest: RITOOLS_TIMESTAMP
	R -q -e "library(RItools, lib.loc = '.local')" \
			 -e "library(SparseM)"  \
		   -e "library(testthat)" \
			 -e "auto_test_package('.')"

build:
	R --vanilla CMD Build .

check: build
	R --vanilla CMD Check RItools_0.1-12.tar.gz

clean:
	git clean

