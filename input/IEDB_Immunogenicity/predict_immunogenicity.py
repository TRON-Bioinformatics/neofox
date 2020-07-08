#!/usr/bin/env python

'''
Created on 08.10.2015
@author: Dorjee Gyaltsen
'''

import os
import sys
from optparse import OptionParser

from logzero import logger

from input.helpers import intermediate_files


class Prediction():

    def isint(self, x):
        try:
            a = float(x)
            b = int(a)
        except ValueError:
            return False
        else:
            return a == b

    def validate(self, options, args):

        allele_dict = {"H-2-Db": "2,5,9", "H-2-Dd": "2,3,5", "H-2-Kb": "2,3,9", "H-2-Kd": "2,5,9", "H-2-Kk": "2,8,9",
                       "H-2-Ld": "2,5,9", "HLA-A0101": "2,3,9", "HLA-A0201": "1,2,9", "HLA-A0202": "1,2,9",
                       "HLA-A0203": "1,2,9", "HLA-A0206": "1,2,9", "HLA-A0211": "1,2,9", "HLA-A0301": "1,2,9",
                       "HLA-A1101": "1,2,9", "HLA-A2301": "2,7,9", "HLA-A2402": "2,7,9", "HLA-A2601": "1,2,9",
                       "HLA-A2902": "2,7,9", "HLA-A3001": "1,3,9", "HLA-A3002": "2,7,9", "HLA-A3101": "1,2,9",
                       "HLA-A3201": "1,2,9", "HLA-A3301": "1,2,9", "HLA-A6801": "1,2,9", "HLA-A6802": "1,2,9",
                       "HLA-A6901": "1,2,9", "HLA-B0702": "1,2,9", "HLA-B0801": "2,5,9", "HLA-B1501": "1,2,9",
                       "HLA-B1502": "1,2,9", "HLA-B1801": "1,2,9", "HLA-B2705": "2,3,9", "HLA-B3501": "1,2,9",
                       "HLA-B3901": "1,2,9", "HLA-B4001": "1,2,9", "HLA-B4002": "1,2,9", "HLA-B4402": "2,3,9",
                       "HLA-B4403": "2,3,9", "HLA-B4501": "1,2,9", "HLA-B4601": "1,2,9", "HLA-B5101": "1,2,9",
                       "HLA-B5301": "1,2,9", "HLA-B5401": "1,2,9", "HLA-B5701": "1,2,9", "HLA-B5801": "1,2,9"}

        sequence_text = open(args[0], "r").read().split()
        for peptide in sequence_text:
            for amino_acid in peptide.strip():
                if not amino_acid.upper() in "ACDEFGHIKLMNPQRSTVWY":
                    print(("Sequence: '%s' contains an invalid character: '%c' at position %d." % (
                    peptide, amino_acid, peptide.find(amino_acid))))
                    sys.exit(1)

        custom_mask = options.custom_mask

        if custom_mask:
            custom_mask_list = list(map(int, custom_mask.split(",")))
            if sum(n < 0 for n in custom_mask_list) > 0:
                print("custom-mask should be greater then zero.")
                sys.exit(1)

            max_length = max(custom_mask_list)
            if not all([len(sequence) >= max_length for sequence in sequence_text]):
                print(("custom length '{}' cannot be greater then the peptide length.".format(max_length)))
                sys.exit(1)

        if custom_mask:
            if not all(self.isint(num) for num in custom_mask.split(",")):
                print("custom-mask should be a single number or comma-separated list of numbers.")
                sys.exit(1)

        allele = options.allele
        allele = allele.replace("*", "").replace(":", "") if allele else None

        # Check if both custom_mask and allele options given
        if custom_mask and allele:
            print((
                      "* Allele {} has default value {}.\n* When both \'custom_mask\' and \'allele\' options are used, latter takes precedence over former.\n".format(
                          allele, allele_dict[allele], '--custom_mask')))

        # Check if allele is included in the available alleles
        if allele in allele_dict:
            custom_mask = allele_dict[allele]

        # Check if allele option is used and is in the available alleles
        if allele:
            if allele not in allele_dict:
                print(("Allele {} is not available.".format(allele)))
                sys.exit(1)

        return sequence_text, custom_mask, allele

    def predict(self, cleaned_data):
        '''Returns the prediction result.'''

        immunoscale = {"A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380, "G": 0.110, "H": 0.105, "I": 0.432,
                       "K": -0.700, "L": -0.036, "M": -0.570, "N": -0.021, "P": -0.036, "Q": -0.376, "R": 0.168,
                       "S": -0.537, "T": 0.126, "V": 0.134, "W": 0.719, "Y": -0.012}
        immunoweight = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]

        result_list = []

        sequence_text, custom_mask, allele = cleaned_data

        for pep in sequence_text:
            peptide = pep.upper()
            peplen = len(peptide)

            cterm = peplen - 1
            score = 0
            count = 0

            if not custom_mask:
                mask_num = [0, 1, cterm]
                mask_out = [1, 2, "cterm"]
            elif custom_mask:
                try:
                    mask_str = custom_mask.split(",")
                    mask_num = list(map(int, mask_str))
                    mask_num = [x - 1 for x in mask_num]
                    mask_out = [x + 1 for x in mask_num]
                except IOError as e:
                    logger.exception(e)
            else:
                self.mask_num = []
                self.mask_out = [1, 2, "cterm"]

            if peplen > 9:
                pepweight = immunoweight[:5] + ((peplen - 9) * [0.30]) + immunoweight[5:]
            else:
                pepweight = immunoweight

            try:
                for pos in peptide:
                    if pos not in list(immunoscale.keys()):
                        raise KeyError()
                    elif count not in mask_num:
                        score += pepweight[count] * immunoscale[pos]
                        count += 1
                    else:
                        count += 1
                result_list.append([peptide, len(peptide), round(score, 5)])

            except IOError as e:
                logger.exception(e)
            except Exception as e:
                logger.exception(e)
                raise

        # Sort by the last column value (score)       
        result_list.sort(key=lambda tup: tup[-1], reverse=True)

        # Column headers for the result
        header_list = ['peptide', 'length', 'score']
        result_list.insert(0, header_list)

        if allele:
            print("allele: {}".format(allele))
        print("masking: {0}".format('custom' if custom_mask else 'default'))
        print("masked variables: {0}\n".format(mask_out))

        # Print comma separate result
        for result in result_list:
            result = list(map(str, result))
            print(','.join(result))

    def create_csv(self, mask_choice, mask_out, data):
        import csv

        tmpdir = './output'

        # Create a temporary file inside the tmp/ directory
        tmpfile = intermediate_files.create_temp_file(prefix="immunogenicity_", suffix=".csv", dir=tmpdir)
        with open(tmpfile, 'wb') as result:
            writer = csv.writer(result, delimiter=',')
            data.insert(0, ['masking: ', '{0}'.format(mask_choice)])
            data.insert(1, ['masked variables: ', '{0}'.format(mask_out)])
            for score in data:
                writer.writerow(score)
        return tmpfile

    def commandline_help(self):
        print("""
* Available options:
--------------------
custom_mask or allele
 
1) To use the tool with default parameters:
-------------------------------------------
python predict_immunogenicity.py [input-file]
Example: python predict_immunogenicity.py example/test.txt
 
2) To use the tool with custom_mask option:
-------------------------------------------
python predict_immunogenicity.py [options] [input-file]
Example: python predict_immunogenicity.py --custom_mask=2,3,9 example/test.txt

3) To use the tool with allele option:
--------------------------------------
python predict_immunogenicity.py [options] [input-file]
Example: python predict_immunogenicity.py --allele=HLA-A01:01 example/test.txt

4) To list all available alleles:
---------------------------------
python predict_immunogenicity.py --allele_list

You can also use help option (-h or --help) for more information:
-----------------------------------------------------------------
python predict_immunogenicity.py --help
""")

    def commandline_allele(self):
        '''Return all available alleles.'''
        print("""
========================================================================= 
| Alleles available for Class-I Immunogenicity:                         | 
|-----------------------------------------------------------------------| 
| H-2-Db    | H-2-Dd    | H-2-Kb    | H-2-Kd    | H-2-Kk    | H-2-Ld    | 
| HLA-A0101 | HLA-A0201 | HLA-A0202 | HLA-A0203 | HLA-A0206 | HLA-A0211 | 
| HLA-A0301 | HLA-A1101 | HLA-A2301 | HLA-A2402 | HLA-A2601 | HLA-A2902 | 
| HLA-A3001 | HLA-A3002 | HLA-A3101 | HLA-A3201 | HLA-A3301 | HLA-A6801 | 
| HLA-A6802 | HLA-A6901 | HLA-B0702 | HLA-B0801 | HLA-B1501 | HLA-B1502 | 
| HLA-B1801 | HLA-B2705 | HLA-B3501 | HLA-B3901 | HLA-B4001 | HLA-B4002 | 
| HLA-B4402 | HLA-B4403 | HLA-B4501 | HLA-B4601 | HLA-B5101 | HLA-B5301 | 
| HLA-B5401 | HLA-B5701 | HLA-B5801                                     |
-------------------------------------------------------------------------
""")

    def main(self):
        parser = OptionParser(usage="usage: %prog [options] input", version="%prog 1.1")

        parser.add_option('-m', '--custom_mask',
                          help="Specify comma separated numbers for additional anchor positions (default: 1st, 2nd, and C-terminus amino acids).")

        parser.add_option('-a', '--allele',
                          help="Specify an allele.")

        parser.add_option('-l', '--allele_list',
                          action='store_true',
                          default=False,
                          help="List all available alleles.")

        (options, args) = parser.parse_args()

        if options.allele_list:
            self.commandline_allele()
            sys.exit(0)

        if len(args) == 0:
            parser.error("Input file is not specified.")

        if os.path.isfile(args[0]) == False:
            self.commandline_help()
            sys.exit(0)

        self.predict(self.validate(options, args))


if __name__ == '__main__':
    pred = Prediction()
    pred.main()
