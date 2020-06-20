import random
import tempfile

import dotenv
from Bio.Alphabet.IUPAC import IUPACData

from input.references import ReferenceFolder, DependenciesConfiguration


def load_references():
    dotenv.load_dotenv()
    return ReferenceFolder(), DependenciesConfiguration()


def create_temp_aminoacid_fasta_file():
    fastafile = tempfile.NamedTemporaryFile(mode='w', delete=False)
    with fastafile as f:
        f.write("".join(random.choices(list(IUPACData.protein_letters), k=25)))
    return fastafile
