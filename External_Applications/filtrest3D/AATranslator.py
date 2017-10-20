# Copyright '2006 Michal Pietal
# See license.txt

from Bio import SeqUtils
from Bio.SCOP import Raf
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein #,IUPACProtein,

class AATranslator:

    """The object for transforming aminoacid three-letter <--> one-letter sequence representation.
    
    get1letter      -   three-letter code to single letter translator.
    get1letter_seq  -   sequence version of a get1letter.
    get1original
    get1original_seq
    get3letter      -   amino acid letter to three letter code translator.
    get3letter_seq  -   sequence version of a get3letter.
    get33original
    get33original_seq
    get_alphabet    -   returns the aminoacid alphabet as a string.
    get_codes       -   returns ordered list of aminoacid codes as a list.
    get_codes_extended
    get_letters     -   returns ordered aminoacid alphabet as a list.
    make_ungapped
    """
    
    def __init__(self, unknown1letter='X', unknown3letter='UNK', gapchar='', separator='', UPPER=1):
        """Optional parameters are: unknown letter and code (when unknown code or letter encountered, respecively), 
        gap code for single letter, separator between codes in returned three-letter sequence, plus UPPER flag whether
        not to distinguish upper/lower case as an input - output is always uppercase.
        """
        self.alphabet = ExtendedIUPACProtein()
        self.unknown1 = unknown1letter 
        self.unknown3 = unknown3letter
        self.gap_char = gapchar
        self.separator = separator
        self.upper=UPPER
        
        self.blankseq3=['Xer','Xaa','Ter','Sel']
    
    def get1letter(self, aa_code):
        """Translation from three-letter code to aminoacid letter.  Faster than get1letter_seq."""
        if self.upper == 1:
            aa_code = aa_code.upper()
        aa = self.unknown1
        for aa_letter in self.alphabet.letters:
            if SeqUtils.seq3(aa_letter).upper() == aa_code:
                aa = aa_letter
                break
        aa = aa.upper()
        return aa
        
    def get1letter_seq(self, aa_code_string):
        """Translation from three-letter string to aminoacid single-letter sequence."""
        sequence = ''
        len_sep = len(self.separator)
        i=0
        while i<len(aa_code_string):
            aa_code = aa_code_string[i:i+3]
            if aa_code == 3*self.gap_char:
                aa = self.gap_char
            else:
                if self.upper == 1:
                    aa_code = aa_code.upper()
                aa = self.unknown1
                for aa_letter in self.alphabet.letters:
                    if SeqUtils.seq3(aa_letter).upper() == aa_code:
                        aa = aa_letter
                        break
            sequence = sequence + aa
            i = i + 3
            if aa_code_string[i:i+len_sep] == self.separator:
                i = i + len_sep
        sequence = sequence.upper()
        return sequence
    
    def get1original(self, aa_code):
        """."""
        if self.upper == 1:
            aa_code = aa_code.upper()
        if Raf.to_one_letter_code.has_key(aa_code):
            aa = Raf.to_one_letter_code[aa_code]
            return aa
        else:
            return self.unknown1
    
    def get1original_seq(self, aa_code_string):
        """."""
        sequence = ''
        len_sep = len(self.separator)
        i=0
        while i<len(aa_code_string):
            aa_code = aa_code_string[i:i+3]
            if aa_code == 3*self.gap_char:
                aa = self.gap_char
            else:
                if self.upper == 1:
                    aa_code = aa_code.upper()
                if Raf.to_one_letter_code.has_key(aa_code):
                    aa = Raf.to_one_letter_code[aa_code]
                else:
                    aa = self.unknown1
            sequence = sequence + aa
            i = i + 3
            if aa_code_string[i:i+len_sep] == self.separator:
                i = i + len_sep
        sequence = sequence.upper()
        return sequence
        
    def get3letter(self, aa):
        """Translation from aminoacid letter to three-letter code.  Faster than get3letter_seq."""
        if len(aa) > 1:
            return self.unknown3
        if self.upper == 1:
            aa = aa.upper()
        code = SeqUtils.seq3(aa)
        if code in self.blankseq3:
            code = self.unknown3
        code = code.upper()
        return code
        
    def get3letter_seq(self, aa_string):
        """Simple translation from aminoacid letter string to three-letter string."""
        if self.upper == 1:
            aa_string = aa_string.upper()
        code = ''
        n=len(aa_string)
        for i in range(n):
            if aa_string[i] == self.gap_char:
                code = code + 3*self.gap_char
            else:
                code_letter = SeqUtils.seq3(aa_string[i])
                if code_letter in self.blankseq3:
                    code = code + self.unknown3
                else:
                    code = code + code_letter
            if i < n-1:
                code = code + self.separator
        code = code.upper()
        return code
    
    def get33original(self, aa_code):
        """."""
        if self.upper == 1:
            aa_code = aa_code.upper()
        if Raf.to_one_letter_code.has_key(aa_code):
            aa_code_original = SeqUtils.seq3(Raf.to_one_letter_code[aa_code])
            if aa_code_original in self.blankseq3:
                return self.unknown3
            else:
                return aa_code_original.upper()
        else:
            return self.unknown3
    
    def get33original_seq(self, aa_code_string):
        """."""
        if self.upper == 1:
            aa_code_string = aa_code_string.upper()
        code_original = ''
        len_sep = len(self.separator)
        i=0
        while i<len(aa_code_string):
            aa_code = aa_code_string[i:i+3]
            if aa_code == 3*self.gap_char:
                aa_code_original = 3*self.gap_char
            elif aa_code == self.unknown3:
                aa_code_original = self.unknown3
            else:
                if Raf.to_one_letter_code.has_key(aa_code):
                    aa_code_original = SeqUtils.seq3(Raf.to_one_letter_code[aa_code])
                    if aa_code_original in self.blankseq3:
                        aa_code_original = self.unknown3
                else:
                    aa_code_original = self.unknown3
            code_original = code_original + aa_code_original
            i = i + 3
            if aa_code_string[i:i+len_sep] == self.separator:
                i = i + len_sep
                code_original = code_original + self.separator
        code_original = code_original.upper()
        return code_original
    
    def get_alphabet(self):
        """Returns the aminoacid alphabet used."""
        return self.alphabet
    
    def get_codes(self):
        """Returns the ordered list of aminoacid codes."""
        c = []
        for letter in self.alphabet.letters:
            if letter == 'X':
                c.append(self.unknown3)
            else:
                c.append(SeqUtils.seq3(letter).upper())
        return c
    
    def get_codes_extended(self):
        """."""
        ce = Raf.to_one_letter_code.keys()
        ce[ce.index('UNK')] = self.unknown3
        return ce
    
    def get_letters(self):
        """Returns the ordered list of aminoacid letters."""
        l = []
        for letter in self.alphabet.letters:
            if letter == 'X':
                l.append(self.unknown1)
            else:
                l.append(letter)
        return l
    
    def make_ungapped(self, sequence):
        """."""
        ungapped_sequence = ''
        for letter in sequence:
            if letter != self.gap_char:
                ungapped_sequence = ungapped_sequence + letter
        return ungapped_sequence

if __name__=="__main__":

    from sys import argv
    
    mode = argv[1]
    arg = argv[2]
    
    if mode == '-3':
        print AATranslator().get1letter_seq(arg)
    elif mode == '-1':
        print AATranslator().get3letter_seq(arg)
