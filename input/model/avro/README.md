
Download avro tools from https://apache.lauf-forum.at/avro/avro-1.9.2/java/avro-tools-1.9.2.jar

Convert the neoepitope.avdl files into neoepitope.avpr:
```
java -jar ~/bin/avro-tools-1.9.2.jar idl neoepitope.avdl neoepitope.avpr
```

Convert the neoepitope.avdl files into Neoepitope.avsc:
```
java -jar ~/bin/avro-tools-1.9.2.jar idl2schemata neoepitope.avdl
```

Installed avrodoc following instructions here https://github.com/ept/avrodoc

Generate Python code using the following script:
```
import json
from avrogen import write_schema_files


schema_json = json.dumps(json.load(open('Neoepitope.avsc', 'r')))
output_directory = '../input/model'
write_schema_files(schema_json, output_directory)
```



