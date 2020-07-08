import json
import os.path
import decimal
import datetime
import six
from avrogen.dict_wrapper import DictWrapper
from avrogen import avrojson
from avro import schema as avro_schema
if six.PY3:    from avro.schema import SchemaFromJSONData as make_avsc_object
    
else:
    from avro.schema import make_avsc_object
    


def __read_file(file_name):
    with open(file_name, "r") as f:
        return f.read()

def __get_names_and_schema(file_name):
    names = avro_schema.Names()
    schema = make_avsc_object(json.loads(__read_file(file_name)), names)
    return names, schema

__NAMES, SCHEMA = __get_names_and_schema(os.path.join(os.path.dirname(__file__), "schema.avsc"))
__SCHEMAS = {}
def get_schema_type(fullname):
    return __SCHEMAS.get(fullname)
__SCHEMAS = dict((n.fullname.lstrip("."), n) for n in six.itervalues(__NAMES.names))


class SchemaClasses(object):
    
    
    pass
    class neoantigen(object):
        
        class GeneClass(DictWrapper):
            
            """
            
            """
            
            
            RECORD_SCHEMA = get_schema_type("neoantigen.Gene")
            
            
            def __init__(self, inner_dict=None):
                super(SchemaClasses.neoantigen.GeneClass, self).__init__(inner_dict)
                if inner_dict is None:
                    self.assembly = SchemaClasses.neoantigen.GeneClass.RECORD_SCHEMA.fields[0].default
                    self.gene = None
                    self.transcriptIdentifier = str()
            
            
            @property
            def assembly(self):
                """
                :rtype: str
                """
                return self._inner_dict.get('assembly')
            
            @assembly.setter
            def assembly(self, value):
                #"""
                #:param str value:
                #"""
                self._inner_dict['assembly'] = value
            
            
            @property
            def gene(self):
                """
                :rtype: str
                """
                return self._inner_dict.get('gene')
            
            @gene.setter
            def gene(self, value):
                #"""
                #:param str value:
                #"""
                self._inner_dict['gene'] = value
            
            
            @property
            def transcriptIdentifier(self):
                """
                :rtype: str
                """
                return self._inner_dict.get('transcriptIdentifier')
            
            @transcriptIdentifier.setter
            def transcriptIdentifier(self, value):
                #"""
                #:param str value:
                #"""
                self._inner_dict['transcriptIdentifier'] = value
            
            
        class MutationClass(DictWrapper):
            
            """
            
            """
            
            
            RECORD_SCHEMA = get_schema_type("neoantigen.Mutation")
            
            
            def __init__(self, inner_dict=None):
                super(SchemaClasses.neoantigen.MutationClass, self).__init__(inner_dict)
                if inner_dict is None:
                    self.position = int()
                    self.wildTypeAminoacid = str()
                    self.mutatedAminoacid = str()
                    self.leftFlankingRegion = str()
                    self.sizeLeftFlankingRegion = None
                    self.rightFlankingRegion = str()
                    self.sizeRightFlankingRegion = None
            
            
            @property
            def position(self):
                """
                :rtype: int
                """
                return self._inner_dict.get('position')
            
            @position.setter
            def position(self, value):
                #"""
                #:param int value:
                #"""
                self._inner_dict['position'] = value
            
            
            @property
            def wildTypeAminoacid(self):
                """
                :rtype: str
                """
                return self._inner_dict.get('wildTypeAminoacid')
            
            @wildTypeAminoacid.setter
            def wildTypeAminoacid(self, value):
                #"""
                #:param str value:
                #"""
                self._inner_dict['wildTypeAminoacid'] = value
            
            
            @property
            def mutatedAminoacid(self):
                """
                :rtype: str
                """
                return self._inner_dict.get('mutatedAminoacid')
            
            @mutatedAminoacid.setter
            def mutatedAminoacid(self, value):
                #"""
                #:param str value:
                #"""
                self._inner_dict['mutatedAminoacid'] = value
            
            
            @property
            def leftFlankingRegion(self):
                """
                :rtype: str
                """
                return self._inner_dict.get('leftFlankingRegion')
            
            @leftFlankingRegion.setter
            def leftFlankingRegion(self, value):
                #"""
                #:param str value:
                #"""
                self._inner_dict['leftFlankingRegion'] = value
            
            
            @property
            def sizeLeftFlankingRegion(self):
                """
                :rtype: int
                """
                return self._inner_dict.get('sizeLeftFlankingRegion')
            
            @sizeLeftFlankingRegion.setter
            def sizeLeftFlankingRegion(self, value):
                #"""
                #:param int value:
                #"""
                self._inner_dict['sizeLeftFlankingRegion'] = value
            
            
            @property
            def rightFlankingRegion(self):
                """
                :rtype: str
                """
                return self._inner_dict.get('rightFlankingRegion')
            
            @rightFlankingRegion.setter
            def rightFlankingRegion(self, value):
                #"""
                #:param str value:
                #"""
                self._inner_dict['rightFlankingRegion'] = value
            
            
            @property
            def sizeRightFlankingRegion(self):
                """
                :rtype: int
                """
                return self._inner_dict.get('sizeRightFlankingRegion')
            
            @sizeRightFlankingRegion.setter
            def sizeRightFlankingRegion(self, value):
                #"""
                #:param int value:
                #"""
                self._inner_dict['sizeRightFlankingRegion'] = value
            
            
        class NeoantigenClass(DictWrapper):
            
            """
            A neoantigen minimal definition
            """
            
            
            RECORD_SCHEMA = get_schema_type("neoantigen.Neoantigen")
            
            
            def __init__(self, inner_dict=None):
                super(SchemaClasses.neoantigen.NeoantigenClass, self).__init__(inner_dict)
                if inner_dict is None:
                    self.gene = SchemaClasses.neoantigen.GeneClass()
                    self.mutation = SchemaClasses.neoantigen.MutationClass()
                    self.expressionValue = None
                    self.clonalityEstimation = None
                    self.variantAlleleFrequency = None
            
            
            @property
            def gene(self):
                """
                :rtype: SchemaClasses.neoantigen.GeneClass
                """
                return self._inner_dict.get('gene')
            
            @gene.setter
            def gene(self, value):
                #"""
                #:param SchemaClasses.neoantigen.GeneClass value:
                #"""
                self._inner_dict['gene'] = value
            
            
            @property
            def mutation(self):
                """
                :rtype: SchemaClasses.neoantigen.MutationClass
                """
                return self._inner_dict.get('mutation')
            
            @mutation.setter
            def mutation(self, value):
                #"""
                #:param SchemaClasses.neoantigen.MutationClass value:
                #"""
                self._inner_dict['mutation'] = value
            
            
            @property
            def expressionValue(self):
                """
                :rtype: float
                """
                return self._inner_dict.get('expressionValue')
            
            @expressionValue.setter
            def expressionValue(self, value):
                #"""
                #:param float value:
                #"""
                self._inner_dict['expressionValue'] = value
            
            
            @property
            def clonalityEstimation(self):
                """
                :rtype: bool
                """
                return self._inner_dict.get('clonalityEstimation')
            
            @clonalityEstimation.setter
            def clonalityEstimation(self, value):
                #"""
                #:param bool value:
                #"""
                self._inner_dict['clonalityEstimation'] = value
            
            
            @property
            def variantAlleleFrequency(self):
                """
                :rtype: float
                """
                return self._inner_dict.get('variantAlleleFrequency')
            
            @variantAlleleFrequency.setter
            def variantAlleleFrequency(self, value):
                #"""
                #:param float value:
                #"""
                self._inner_dict['variantAlleleFrequency'] = value
            
            
        pass
        
__SCHEMA_TYPES = {
'neoantigen.Gene': SchemaClasses.neoantigen.GeneClass,
    'neoantigen.Mutation': SchemaClasses.neoantigen.MutationClass,
    'neoantigen.Neoantigen': SchemaClasses.neoantigen.NeoantigenClass,
    
}
_json_converter = avrojson.AvroJsonConverter(use_logical_types=False, schema_types=__SCHEMA_TYPES)

