from unittest import TestCase

from neofox.model.wrappers import AnnotationFactory


class TestAnnotationFactory(TestCase):

    def test_float_values_precision_limit(self):
        self.assertEqual("0.12346", AnnotationFactory.build_annotation("test", 0.123456789).value)
        self.assertEqual("0.12346", AnnotationFactory.build_annotation("test", 0.1234567891234).value)
        self.assertEqual("0.12345", AnnotationFactory.build_annotation("test", 0.12345).value)
        self.assertEqual("0.123", AnnotationFactory.build_annotation("test", 0.123).value)

    def test_none(self):
        self.assertEqual("NA", AnnotationFactory.build_annotation("test", None).value)

    def test_booleans(self):
        self.assertEqual("1", AnnotationFactory.build_annotation("test", True).value)
        self.assertEqual("0", AnnotationFactory.build_annotation("test", False).value)

    def test_integers(self):
        self.assertEqual("123", AnnotationFactory.build_annotation("test", 123).value)
        self.assertEqual("-123", AnnotationFactory.build_annotation("test", -123).value)

    def test_strings(self):
        # TODO: I am wondering if we should transform this to NA. Franziska what do you think?
        self.assertEqual("", AnnotationFactory.build_annotation("test", "").value)
        self.assertEqual("blabla", AnnotationFactory.build_annotation("test", "blabla").value)
