.PHONY: clean clean-test clean-pyc clean-build docs help env
.DEFAULT_GOAL := help
define BROWSER_PYSCRIPT
import os, webbrowser, sys
try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT
BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## Remove all artefacts.

clean-build: ## Remove build artefacts.
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## Remove Python file artefacts.
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## Remove test and coverage artefacts.
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

env: ## Create development environment.
	conda env create --file environment.yml
	pip install -e .

lint: ## Check style with pylint.
	pylint src tests

test: ## Run tests quickly with the default Python
	py.test

coverage: ## Check code coverage quickly with the default Python.
	coverage run --source pyim py.test

		coverage report -m
		coverage html
		$(BROWSER) htmlcov/index.html

docs: ## Generate Sphinx HTML documentation, including API docs.
	rm -rf docs/_build
	sphinx-autobuild docs docs/_build

release: clean ## Package and upload a release.
	python setup.py sdist upload
	python setup.py bdist_wheel upload

dist: clean ## Builds source and wheel package.
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

install: clean ## Install the package to the active Python's site-packages.
	python setup.py install

gh-pages: ## Deploy documentation on github-pages.
	git checkout gh-pages
	find ./* -not -path '*/\.*' -prune -exec rm -r "{}" \;
	git checkout develop docs Makefile src AUTHORS.rst CONTRIBUTING.rst HISTORY.rst README.rst
	git reset HEAD
	(cd docs && make html)
	mv -fv docs/_build/html/* ./
	rm -rf docs Makefile src AUTHORS.rst CONTRIBUTING.rst HISTORY.rst README.rst
	touch .nojekyll
	git add -A
	git commit -m "Generated gh-pages for `git log develop -1 --pretty=short --abbrev-commit`" && git push origin gh-pages ; git checkout develop
