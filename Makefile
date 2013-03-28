R = R_LIBS=.local R --vanilla

R: .local/RItools/INSTALLED
	$(R) -q --no-save 

### Package release scripts ###

VERSION=0.1-12
RELEASE_DATE=`date +%Y-%m-%d`
PKG=RItools_$(VERSION)

# a useful helper for scripts who need to know what the package name is going to be
# use: R CMD INSTALL path/to/RItools/$(cd path/to/RItools && make current)
current: 
	@echo $(PKG).tar.gz

# we depend on the makefile so that updates to the version number will force a rebuild
$(PKG): Makefile R/* tests/* inst/tests/* man/* .Rinstignore
	rm -rf $(PKG)
	rsync -a --exclude-from=.gitignore --exclude=.git* --exclude Makefile \
		--exclude=DESCRIPTION.template --exclude=NAMESPACE.static \
		--exclude=lexicon.txt --exclude=README.md --exclude=checkspelling.R \
		--exclude=RItools.Rcheck \
		. $(PKG)

$(PKG)/DESCRIPTION: $(PKG) DESCRIPTION.template 
	sed s/VERSION/$(VERSION)/ DESCRIPTION.template | sed s/DATE/$(RELEASE_DATE)/ > $(PKG)/DESCRIPTION

$(PKG)/NAMESPACE: $(PKG) $(PKG)/DESCRIPTION NAMESPACE.static .local/roxygen2/INSTALLED
	mkdir -p $(PKG)/man
	$(R) -e "library(roxygen2); roxygenize('$(PKG)')"
	cat NAMESPACE.static >> $(PKG)/NAMESPACE

$(PKG).tar.gz: $(PKG) $(PKG)/DESCRIPTION $(PKG)/NAMESPACE NEWS R/* data/* inst/* man/* tests/*
	$(R) --vanilla CMD build $(PKG)

package: $(PKG).tar.gz

# the spell task doesn't need the tar.gz particularly, but it does need DESCRIPTION and roxygen
spell: package 
	$(R) -q --no-save -e "source('checkspelling.R') ; check_spelling('$(PKG)')"

lexicon.txt: package
	$(R) -q --no-save -e "source('checkspelling.R') ; make_dictionary('$(PKG)')"

# For reasons unknown, the check process has a fit if using just the environment variables R_LIBS 
# so, we also use the -l flag. 
check: $(PKG).tar.gz .local/xtable/INSTALLED .local/SparseM/INSTALLED .local/optmatch/INSTALLED
	$(R) CMD check --as-cran --no-multiarch -l .local $(PKG).tar.gz

release: check spell
	git tag -a $(VERSION)
	@echo "Upload $(PKG).tar.gz to cran.r-project.org/incoming"
	@echo "Email to CRAN@R-project.org, subject: 'CRAN submission RItools $(VERSION)'"


# additional dependencies from CRAN
installpkg = mkdir -p .local ; $(R) -e "install.packages('$(1)', repos = 'http://streaming.stat.iastate.edu/CRAN/')" ; date > .local/$(1)/INSTALLED
	
.local/roxygen2/INSTALLED:
	$(call installpkg,roxygen2)

.local/testthat/INSTALLED:
	$(call installpkg,testthat)

.local/SparseM/INSTALLED:
	$(call installpkg,SparseM)

.local/optmatch/INSTALLED:
	$(call installpkg,optmatch)

.local/xtable/INSTALLED:
	$(call installpkg,xtable)
	
# depend on this file to decide if we need to install the local version
.local/RItools/INSTALLED: $(PKG).tar.gz .local/SparseM/INSTALLED .local/optmatch/INSTALLED .local/xtable/INSTALLED
	mkdir -p .local
	$(R) CMD INSTALL --no-multiarch --library=.local $(PKG).tar.gz
	echo `date` > .local/RItools/INSTALLED

# If we run into problems documenting S4 code, this helped
# .local/roxygen2/INSTALLED:
# 	mkdir -p .local
# 	R_LIBS=.local R --vanilla -e "library(devtools) ; install_github(repo = 'roxygen', user = 'klutometis', branch = 's4',args=c('--no-multiarch'))"
# 	echo `date` > .local/roxygen2/INSTALLED

# test is just the internal tests, not the full R CMD Check
test: .local/RItools/INSTALLED .local/testthat/INSTALLED
	$(R) -e "library(RItools, lib.loc = '.local'); library(testthat); test_package('RItools')"

# removes local files, leaves external libraries
clean:
	mv .local .local-clean
	git clean -Xfd
	mv .local-clean .local

clean-deps:
	rm -rf .local
