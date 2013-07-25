R = R_LIBS=.local R

R: .local/RItools/INSTALLED
	$(R) -q --no-save 

### Package release scripts ###

VERSION=0.1-12
RELEASE_DATE=`date +%Y-%m-%d`
PKG=RItools_$(VERSION)

# we depend on the makefile so that updates to the version number will force a rebuild
$(PKG): Makefile R/* tests/* inst/tests/* man/* 
	rm -rf $(PKG)
	rsync -a --exclude-from=.gitignore --exclude=.git* --exclude Makefile \
		--exclude=DESCRIPTION.template --exclude=NAMESPACE.static \
		--exclude=lexicon.txt --exclude=README.md --exclude=checkspelling.R \
		--exclude=RItools.Rcheck \
		. $(PKG)

$(PKG)/DESCRIPTION: $(PKG) DESCRIPTION.template 
	sed s/VERSION/$(VERSION)/ DESCRIPTION.template | sed s/DATE/$(RELEASE_DATE)/ > $(PKG)/DESCRIPTION

$(PKG)/NAMESPACE: $(PKG) $(PKG)/DESCRIPTION NAMESPACE.static 
	mkdir -p $(PKG)/man
	$(R) -e "library(roxygen2); roxygenize('$(PKG)')"
	cat NAMESPACE.static >> $(PKG)/NAMESPACE

$(PKG).tar.gz: $(PKG) $(PKG)/DESCRIPTION $(PKG)/NAMESPACE NEWS R/* data/* inst/* man/* tests/*
	R --vanilla CMD build $(PKG)

package: $(PKG).tar.gz

# the spell task doesn't need the tar.gz particularly, but it does need DESCRIPTION and roxygen
spell: package 
	$(R) -q --no-save -e "source('checkspelling.R') ; check_spelling('$(PKG)')"

lexicon.txt: package
	$(R) -q --no-save -e "source('checkspelling.R') ; make_dictionary('$(PKG)')"

# we don't use $(R) here in case the locally installed libs interfere with the check process
check: $(PKG).tar.gz 
	R --vanilla CMD check --as-cran --no-multiarch $(PKG).tar.gz

release: check spell
	git tag -a $(VERSION)
	@echo "Upload $(PKG).tar.gz to cran.r-project.org/incoming"
	@echo "Email to CRAN@R-project.org, subject: 'CRAN submission RItools $(VERSION)'"

# depend on this file to decide if we need to install the local version
.local/RItools/INSTALLED: $(PKG).tar.gz
	mkdir -p .local
	R --vanilla CMD INSTALL --no-multiarch --library=.local $(PKG).tar.gz
	echo `date` > .local/RItools/INSTALLED

# If we run into problems documenting S4 code, this helped
# .local/roxygen2/INSTALLED:
# 	mkdir -p .local
# 	R_LIBS=.local R --vanilla -e "library(devtools) ; install_github(repo = 'roxygen', user = 'klutometis', branch = 's4',args=c('--no-multiarch'))"
# 	echo `date` > .local/roxygen2/INSTALLED

# test is just the internal tests, not the full R CMD Check
test: .local/RItools/INSTALLED
	R_LIBS=.local R --vanilla -q -e "library(RItools, lib.loc = '.local'); library(testthat); test_package('RItools')"

clean:
	git clean -xfd

