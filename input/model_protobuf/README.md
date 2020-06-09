
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




