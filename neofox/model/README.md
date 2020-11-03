# Data models

## Compiling the models and generating source code

Install the protocol buffers compiler. Download and unzip from the releases page https://github.com/protocolbuffers/protobuf/releases/tag/v3.12.3

Install the tool betterproto to make generated Python source code a bit better:
```
pip install "betterproto[compiler]"
```

Generate Python source code from neoepitope.proto:
```
protoc -I . --python_betterproto_out=. neoantigen.proto
```

This will generate `neoantigen.py` and `__init__.py` files.


## Documentation

Documentation is generated using this tool https://github.com/pseudomuto/protoc-gen-doc.

Run the following to generate markdown documentation:
```
sudo docker run --rm   -v $(pwd):/out   -v $(pwd):/protos pseudomuto/protoc-gen-doc --doc_opt=/protos/models_template.tmpl,models.md
```

Using `/protos/models_template.tmpl` template allows us to customise the markdown output.

To integrate this documentation into the package documentation copy it into the right folder:
```
cp models.md ../../docs/source/07_models.md
```


