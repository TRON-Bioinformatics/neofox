import dotenv
import tempfile
from input.references import ReferenceFolder
from Bio.Alphabet.IUPAC import IUPACData
import random


def load_references():
    dotenv.load_dotenv()
    return ReferenceFolder()


def create_temp_aminoacid_fasta_file():
    fastafile = tempfile.NamedTemporaryFile(mode='w', delete=False)
    with fastafile as f:
        # TODO: change this to random.choices(k=25) once in Python 3
        f.write("".join([random.choice(list(IUPACData.protein_letters)) for _ in range(25)]))
    return fastafile