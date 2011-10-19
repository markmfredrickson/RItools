# the default task starts R, with RItools loaded (see .Rprofile)
interactive: .local/RItools
	R --quiet --no-save   

.local:
	mkdir .local

.local/RItools: R/* tests/* .local
	R --vanilla CMD Install --library=.local .

autotest: .local/RItools/
	R -q -e "library(RItools, lib.loc = '.local')" \
			 -e "library(SparseM)"  \
		   -e "library(testthat)" \
			 -e "auto_test_package('.')"


