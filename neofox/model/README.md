

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

Documentation is generated using this tool https://github.com/pseudomuto/protoc-gen-doc. Follow the instructions in the README for installation and then add the binary to your PATH.

Run:
```
protoc --doc_out=. --doc_opt=html,../../docs/build/html/neoantigen.html neoantigen.proto
```


