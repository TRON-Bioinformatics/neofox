# Data models

## Compiling the models and generating source code

Install the protocol buffers compiler. Download and unzip from the releases page https://github.com/protocolbuffers/protobuf/releases/tag/v3.12.3

Install the tool betterproto to make generated Python source code a bit better:
```
pip install "betterproto[compiler]"
```

Generate Python source code from neoepitope.proto:
```
make models
```

This will generate `neoantigen.py` and `__init__.py` files.


## Documentation

Documentation is generated using this tool https://github.com/pseudomuto/protoc-gen-doc through a docker instance, but this requires sudo permissions.

Run the following to generate markdown documentation:
```
make html
```

Using `/protos/models_template.tmpl` template allows us to customise the markdown output.
