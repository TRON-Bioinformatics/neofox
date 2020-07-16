from logzero import logger

immunoweight = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]

immunoscale = {"A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380, "G": 0.110, "H": 0.105, "I": 0.432,
                        "K": -0.700, "L": -0.036, "M": -0.570, "N": -0.021, "P": -0.036, "Q": -0.376, "R": 0.168,
                        "S": -0.537, "T": 0.126, "V": 0.134, "W": 0.719, "Y": -0.012}

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


class IEDBimmunogenicity:

    def predict_immunogenicity(self, pep, allele):

        custom_mask = allele_dict.get(allele, False)
        peptide = pep.upper()
        peplen = len(peptide)

        cterm = peplen - 1
        score = 0
        count = 0

        if not custom_mask:
            mask_num = [0, 1, cterm]
        elif custom_mask:
            mask_str = custom_mask.split(",")
            mask_num = list(map(int, mask_str))
            mask_num = [x - 1 for x in mask_num]

        if peplen > 9:
            pepweight = immunoweight[:5] + ((peplen - 9) * [0.30]) + immunoweight[5:]
        else:
            pepweight = immunoweight
        try:
            for pos in peptide:
                if pos not in list(immunoscale.keys()):
                    logger.debug("{} {} {}".format(pos, pep, allele))
                    raise KeyError()
                elif count not in mask_num:
                    score += pepweight[count] * immunoscale[pos]
                    count += 1
                else:
                    count += 1
        except Exception as ex:
            logger.exception(ex)
            raise ex
        return score

    def calc_IEDB_immunogenicity(self, epitope, mhc_allele, mhc_score, affin_filtering=False):
        '''
        This function determines the IEDB immunogenicity score
        '''
        score = "NA"
        try:
            if affin_filtering:
                if float(mhc_score) < 500.0:
                    score = self.predict_immunogenicity(
                        epitope, mhc_allele.replace("*", "").replace(":",
                                                                     ""))
            else:
                score = self.predict_immunogenicity(
                    epitope, mhc_allele.replace("*", "").replace(":",
                                                                 ""))
        except ValueError:
            score = "NA"
        return score


