#R Package development system
#By Andrew G Hart (2014)


#Include package-specific variables, etc.
include ThisPackage

#Get the package version from the DESCRIPTION file and determine the name of the package tar ball
PKG-VER=$(shell cat $(PKG-NAME)/DESCRIPTION |grep "Version:" |sed 's/Version: //')
PKG-TARBALL=$(PKG-NAME)_$(PKG-VER).tar.gz

#Determine files that belong to the package.
#If any of the files in SRC-FILES change, then make or make build will cause the tar ball to be rebuilt.
basefiles=DESCRIPTION NAMESPACE INDEX configure cleanup LICENSE LICENCE \
NEWS INSTALL README README.md ChangeLog configure.win cleanup.win BinaryFiles \
R/* R/unix/* R/windows/* man/*.Rd man/unix/*.Rd man/windows/*.Rd \
data/* demo/* exec/* inst/* po/* tests/* tools/* vignettes/* src/*
allfiles=$(foreach f,$(basefiles),$(wildcard $(PKG-NAME)/$(f)))
excludefiles=$(PKG-NAME)/$(PKG-NAME)-Ex.R $(EXCLUDE-FILES)
SRC-FILES=$(filter-out $(excludefiles),$(allfiles))

#Some notes about allowable files in package directory structure
#R: *.[qrsRS], *.in sysdata.rda, unix/*.[qrsRS], windows/*.[qrsRS]
#man: *.[rR]d, unix/*.[rR]d, windows/*.[rR]d
#src: Makevars, Makevars.win, Makefile, Makefile.win, *.c, *.cc, *.cpp, *.h, *.f, *.f90, *.f95, *.m, *.mm, install.libs.R
#demo: *.r, *.R, 00Index
#tests: *.R, *.Rin, *.Rout.save, Examples/pkg-Ex.Rout.save

#Not necessary, but define the args variable as empty initially
args=

#Determine operating system
ifdef OS
  ifneq ($(findstring Windows,$(OS)),)
    OSType=Windows
  else
    OSType=$(OS)
  endif
else
  OSType=Unix
endif

#Get the name of the pager
ifeq ($OSType),Windows)
  PAGER=$(WINDOWS-PAGER)
else
  ifeq ($(OSType),Unix)
    PAGER=$(UNIX-PAGER)
  else
PAGER=$(WINDOWS-PAGER)
  endif
endif

#Targets

.PHONEY: build
build: $(PKG-TARBALL)

$(PKG-TARBALL): $(SRC-FILES)
	@R CMD build $(PKG-NAME)

.PHONEY: check
check: clean
	@R CMD check $(PKG-NAME) $(args)

.PHONEY: install
install: $(PKG-TARBALL)
	@R CMD INSTALL $(PKG-TARBALL) $(args)

.PHONEY: remove
remove:
	@R CMD REMOVE $(PKG-NAME)
	
.PHONEY: clean
clean:
	@$(RM) -f $(PKG-NAME)/src/*.o
	@$(RM) -f $(PKG-NAME)/src/*.dll
	@$(RM) -f $(PKG-NAME)/src/*.so
	@$(RM) -f $(PKG-NAME)/src/*.rds
	@$(RM) -rf $(PKG-NAME)/src-i386
	@$(RM) -rf $(PKG-NAME)/src-x64
	@$(RM) -f $(PKG-NAME)/$(PKG-NAME)-Ex.R
	@$(RM) -rf $(PKG-NAME).Rcheck
	@find . -name "*.bak" -delete
		
.PHONEY: log
log:
	@$(PAGER) $(PKG-NAME).Rcheck/00check.log

.PHONEY: out
out:
	@$(PAGER) $(PKG-NAME).Rcheck/00install.out

#The echo target is for debugging the PAGER selection and source file code
.PHONEY: echo
echo:
	@echo $(OSType) ... $(PAGER)
	@echo $(SRC-FILES)
	