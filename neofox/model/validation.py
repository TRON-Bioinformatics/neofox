import betterproto

from neofox.model.neoantigen import Neoantigen


class ModelValidator(object):

    @staticmethod
    def validate(model: betterproto.Message):
        # TODO: make this method capture appropriately validation issues whend dealing with int and float
        return model.__bytes__()

    @staticmethod
    def validate_neoantigen(neoantigen: Neoantigen):
        # TODO: check that all fields are consistent
        pass