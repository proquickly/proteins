#!/usr/bin/make -f
# add gitgub action

PROJECT_ROOT="/Users/andy/ws"
PYTHON_VER="3.11.4"

project?=$(notdir $(CURDIR))
PROJECT=$(project)
FULL_PROJECT_PATH = $(PROJECT_ROOT)/$(PROJECT)
visibility?="public"
VISIBILITY=$(visibility)

local: check create initialize addlocalpy installdeps createrepo
docker: create build createrepo

.DEFAULT_GOAL = check 
.PHONY : check create initialize addlocalpy installdeps createrepo scan run clean destroy install test drun help

check:
ifeq ($(PROJECT), pybuild)
	@echo "Must use a parameter when in pybuild dir"
	@echo
endif

help: check
	@echo "---------------------HELP-----------------------"
	@echo "To test the project type make test"
	@echo "To run the project type make run"
	@echo "To update dependencies type make updeps"
	@echo "make project=dirname visibility=puborpri local"
	@echo "To setup a local project enter local"
	@echo "------------------------------------------------"
	@exit 99


create:
	./Makescripts/create $(PROJECT_ROOT) $(PROJECT)

initialize:
	./Makescripts/initialize $(PROJECT_ROOT) $(PROJECT)

addlocalpy:
	./Makescripts/addlocalpy $(PROJECT_ROOT) $(PROJECT) $(PYTHON_VER)

installdeps:
	./Makescripts/installdeps $(PROJECT_ROOT) $(PROJECT)

createrepo:
	./Makescripts/createrepo $(PROJECT_ROOT) $(PROJECT) $(VISIBILITY)

destroy:
	./Makescripts/destroy $(PROJECT_ROOT) $(PROJECT) 

install:
	./install

updeps:
	python -m pip install --upgrade -r requirements.txt 1>/dev/null
	python -m pip install --upgrade -r requirements.dev.txt 1>/dev/null

scan:
	./Makescripts/scan

clean:
	find . -type f -name ‘.DS_Store’ -delete
	find . -type f -name ‘*.pyc’ -delete
	rm -rf __pycache__

run:
	python src/$(PROJECT)/app.py

test:
	pytest -vv tests/

drun:
	doppler run -p $(PROJECT) -c dev -- python src/$(PROJECT)/app.py

pypideploy:
	#

ghadeploy:
	#
