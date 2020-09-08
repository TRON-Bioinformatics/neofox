from neofox.model.neoantigen import Annotation

NOT_AVAILABLE_VALUE = "NA"


class AnnotationFactory(object):

    @staticmethod
    def build_annotation(name, value):
        if isinstance(value, bool):
            # prints booleans as 0/1 strings
            value = "1" if value else "0"
        if not isinstance(value, str) and not isinstance(value, type(None)):
            if isinstance(value, float):
                value = "{0:.5g}".format(round(value, 5))
            else:
                value = str(value)
        if value is None:
            value = NOT_AVAILABLE_VALUE
        return Annotation(name=name, value=value)
