import dotenv
from input.references import ReferenceFolder


def load_references():
    dotenv.load_dotenv()
    return ReferenceFolder()