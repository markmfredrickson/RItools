.local:
	mkdir .local

.local/RItools/: R/* tests/* .local
	R CMD Install --library=.local .

autotest: .local/RItools
	R -q -e "library(RItools, lib.loc = '.local')" \
			 -e "library(testthat)" \
			 -e "auto_test_package('.')"

