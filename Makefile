PROJECT_NAME = ctfr

SOURCE_BASE_FOLDER = src
SOURCE_LOCATION = $(SOURCE_BASE_FOLDER)/$(PROJECT_NAME)
CY_IMPLEMENTATIONS_LOCATIONS = $(SOURCE_LOCATION)/implementations

RM = rm -f
RMR = $(RM) -r
PYTHON_EXEC ?= python3
PIP = python3 -m pip
BUILD = $(PYTHON_EXEC) -m build
SETUP = $(PYTHON_EXEC) setup.py
DIST = dist

TWINE = python3 -m twine
WHEELHOUSE = wheelhouse

install:
	$(PIP) install .

dev:
	$(PIP) install --editable .[dev]

ext:
	$(SETUP) build_ext --inplace

annotate:
	ANNOTATE=1 $(SETUP) build_ext --inplace

uninstall:
	$(PIP) uninstall $(PROJECT_NAME)


