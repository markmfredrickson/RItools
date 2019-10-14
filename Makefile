# load.R fixes a bug with devtool's `help` to enable `help` on
# functions in this package, as well as loading the package
LOAD=R_PROFILE=load.R
RCMD=R -q

interactive:
	@$(LOAD) $(RCMD) --no-save

interactive-emacs:
	@$(LOAD) R_LIBS=.local emacs -nw -f R

.devtools:
	@$(RCMD) -e "devtools:::$(FUNC)($(DEVTOOLSARG))"

DEVTOOLSARG=
dependencies: FUNC=install_deps
dependencies: DEVTOOLSARG=dependencies=TRUE
test: FUNC=test
check: FUNC=check
document: FUNC=document
vignette: FUNC=build_vignettes # To be renabled if we add vignettes
clean-vignette: FUNC=clean_vignettes
build: FUNC=build
dependencies test check document build: .devtools
vignette clean-vignette: .devtools

clean: #clean-vignette
	git clean -Xfd

inst/aspirin_use.rda: nonbuild/aspirin/analyzing_aspirin_data.R \
	nonbuild/aspirin/BPX_H.XPT \
	nonbuild/aspirin/DEMO_H.XPT \
	nonbuild/aspirin/RXQASA_H.XPT
	cd nonbuild/aspirin/ && $(RCMD) -f analyzing_aspirin_data.R && mv analyzing_aspirin_data.rda ../../inst/aspirin_use.rda
