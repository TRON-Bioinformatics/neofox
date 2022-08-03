# Examples

## Demo dataset

Sample input:
- :download:`Neoantigens input file </_static/test_data.tsv>`
- :download:`Patients input file </_static/test_patients.tsv>`

Command:
```
neofox --input-file test_data.tsv --patient-data test_patients.tsv --output-folder your_folder --output-prefix test --with-all-neoepitopes --organism human
```

Sample output:
- :download:`Annotated neoantigens in tabular format </_static/test_neoantigen_candidates_annotated.tsv>`
- :download:`Annotated neoantigens and neoepitopes in JSON format </_static/test_neoantigen_candidates_annotated.json>`
- :download:`Annotated MHC-I neoepitopes in tabular format </_static/test_mhcI_epitope_candidates_annotated.tsv>`
- :download:`Annotated MHC-II neoepitopes in tabular format </_static/test_mhcII_epitope_candidates_annotated.tsv>`