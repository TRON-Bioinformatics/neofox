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
        f.write("".join(random.choices(list(IUPACData.protein_letters), k=25)))
    return fastafile