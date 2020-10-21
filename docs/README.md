
To generate the Python source code automatic documentation run:
```
sphinx-apidoc -f -o docs/source neofox/ neofox/tests/* neofox/references/* neofox/MHC_predictors/* neofox/helpers/* neofox/command_line.py neofox/annotation_resources/* neofox/published_features/* --ext-autodoc --ext-doctest --ext-todo --ext-viewcode
```

This will only generate code documentation for the main neofox API, the models and the exceptions.

To generate the HTML documentation run `make html` from the `docs` folder.