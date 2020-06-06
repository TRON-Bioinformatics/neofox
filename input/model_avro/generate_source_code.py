#!/bin/python
import json
from avrogen import write_schema_files


schema_json = json.dumps(json.load(open('Neoantigen.avsc', 'r')))
output_directory = '.'
write_schema_files(schema_json, output_directory)
