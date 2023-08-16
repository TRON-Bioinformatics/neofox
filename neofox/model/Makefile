models:
	protoc -I . --python_betterproto_out=. neoantigen.proto

html:
	docker run --rm   -v `pwd`:/out   -v `pwd`:/protos pseudomuto/protoc-gen-doc --doc_opt=/protos/models_template.tmpl,models.md
	cp models.md ../../docs/source/05_models.md
