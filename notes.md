setuptools - for creating setup.py
- mainly provides the `setup` function

twine - for distribution, pushing to PyPI
 - primary means of publishing a library

```
# check the status of the distribution
twine check dist/*

# test upload the library
twine upload --repository-url https://test.pypi.org/legacy/ dist/*

# true upload to pypi
twine upload dist/*
```

wheel - for creating the distributable
- uses the `setup.py` file to the create distributable

```
# binary distribution
python setup.py bdist_wheel

# source distribution
python setup.py sdist
```

