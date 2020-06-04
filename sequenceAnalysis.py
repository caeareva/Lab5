#!/usr/bin/env python3

#####################################################################################################
# File: sequenceAnalysis.py 
# Executable: called by a nother class or main function
# Pupose: Program calculates the composition statistics of a genome such as GC content, codon compositions,
#         and codon usage. Program must be executed in std in, otherwise, include the file path in Main().
#
# Student: Carlos Arevalo(caeareva)
# Group: None 
#
#####################################################################################################

'''  
This program calculates properties and statistics of a genome. For example, GC content, codon composition, 
and codon frequency in a genome fasta file. More propertites might be calculate by calling
more classes in sequence analysis. 
'''
import sys
class NucParams:
    '''
    This class uses disctionaries to study and calculate the properties of a fasta file. The class calculates 
    the GC content, amino acid composition and codon usage in a genome. The class uses dictionaries to read and define every three 
    nucleotides as a codon. Then, it finds in the input genome how many times that codon repeats and canculates
    the codon composition and usage.
    '''
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
    }
   
    # creates a dictionary for amino acids, including stop codons
    aaDictionary = {'A' : 0, 'R' : 0, 'N' : 0, 'D' : 0, 'C' : 0, 'E' : 0, 'Q'
                   : 0, 'G' : 0, 'H' : 0, 'I' : 0, 'L' : 0, 'K' : 0, 'M' : 0,
                   'F' : 0, 'P' : 0, 'S' : 0, 'T' : 0, 'W' : 0, 'Y' : 0, 'V' : 0,
                    '-' : 0}
    
    # convert RNA table to DNA table by replacing 'U' by 'T' 
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    
    # creates a dictionary for allowed nucleotides 
    allowedNucleotides = {'A', 'T', 'C', 'G', 'U', 'N'}
    

    def __init__ (self):
        '''
        Function initializes program creating dictionaries, assigning values to codon 
        and nucleotide dictionaries, and stating nucleotide parameters.
        '''
        self.aaComp = NucParams.aaDictionary # access allowedAminos dictioary 
        # method assings a value of zero to values in rnaCodonTable values
        self.codonComp = {key:0 for key in self.rnaCodonTable.keys()}
        # methods assings a value of zero to values in allowedNucleotides values
        self.nucleotideComp = {key:0 for key in self.allowedNucleotides}

    def addSequence (self, inSequence):
        '''
        Function creates an upper case string, translate DNA sequence to RNA sequences,
        and assigns every three nucleotides as a codon adding them to dictionary
        '''
        # If input sequence has lower case nucleotides, method turns them to upper case
        sequence = ''.join(inSequence.split()).upper()
        #sequence = self.cleanSeqCleaner(inSequence)
        for nucleotide in self.allowedNucleotides:
            # counts nucleotides in input sequence
            self.nucleotideComp[nucleotide] += sequence.count(nucleotide)
            pass
        
        # if input file is a DNA sequence, method translates it to RNA
        newSequence = sequence.replace('T', 'U') 
        # iterates and defines codons as three nucleotides in input sequence
        for i in range(0, len(newSequence), 3):
            codons = newSequence[i:i + 3]  
            # if codons in dictionary, add one
            if codons in self.codonComp:
                self.codonComp[codons] += 1
                # defines aa as codons
                aa = self.rnaCodonTable[codons] 
                self.aaComp[aa] += 1 
                pass
            
    def aaComposition(self):
        '''
        Funtion returns self.aaComp dictionary
        '''
        return self.aaComp

    def nucComposition(self):
        '''
        Funtion returns self.nucleotideComp dictionary
        '''
        return self.nucleotideComp

    def codonComposition(self): 
        '''
        Funtion returns self.codonComp dictionary
        '''
        return self.codonComp
    
    def nucCount(self):
        '''
        Function returns the sum of self.nucleotidesComp values
        '''
        return sum(self.nucleotideComp.values())

#####################################################################################################
# Class ProteinParam
#####################################################################################################

"""
This program calculates the the physical and chemical properties of a protein sequence (amino acids).
Inputs has to be a protein sequence, and the program will output the number of amino acids, molecular weights,
molar extinction coefficient, mass extinction coefficients, theoretical isoelectric point (PI), and amino acids composition.
If a character of amino acid is not identified the program then it won't be executable.

Example:
    Protein sequence: VLSPADKTNVKAAW

Output:
    Protein sequence? VLSPADKTNVKAAW
    Number of Amino Acids: 14
    Molecular Weight: 1499.7
    molar Extinction coefficient: 5500.00
    mass Extinction coefficient: 3.67
    Theoretical pI: 9.88
    Amino acid composition:
        A = 21.43%
        C = 0.00%
        D = 7.14%
        E = 0.00%
        F = 0.00%
        G = 0.00%
        H = 0.00%
        I = 0.00%
        K = 14.29%
        L = 7.14%
        M = 0.00%
        N = 7.14%
        P = 7.14%
        Q = 0.00%
        R = 0.00%
        S = 7.14%
        T = 7.14%
        V = 14.29%
        W = 7.14%
        Y = 0.00%
    Protein sequence?
"""

class ProteinParam:
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like
    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    # amino acids molecular weights
    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015  # molecular weight of water
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}  # amino acid absorbances at 280 nm
    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}  # pKa of positively charged Amino Acids
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}  # pKa of negatively charged Amino Acids
    aaNterm = 9.69  # pKa of the Nitrous terminus
    aaCterm = 2.34  # pKa of the Carboxyl terminus

    def __init__(self, protein):
        '''
        Function uses a loop to check allowed characters "Aminos" and then join then
        to create and print a protein string.
        '''
        splitAmino = []  # creates an empty list
        aminos = self.aa2mw.keys()  # access keys in dictionary
        for amino in protein:
            if amino in aminos:
                splitAmino.append(amino)  # appends chars in protein to list
                pass

        # stores and joins aminos seq in object
        newAmino = ''.join(splitAmino).split()
        # stores, converts to upper case, and prints protein string
        self.proteinString = ''.join(newAmino).upper()
        pass

    def aaCount(self):
        '''
        Returns counts of the protein string: amino acids
        '''
        return len(self.proteinString)

    def calculatePI(self):
        '''
        Calculates and return the isoelectric point (PI) of protein sequence
        '''
        lowestCharge = (2**5, 0.00)  # charge in method is an arbitary number
        for ph in range(0, 1400 + 1):  # assigns ph range value from 1 to 14
            pH = ph / 100.  # calculates pH as a float
            charge = self._charge_(pH)
            # if charge is in the range returns lowestcharge
            if charge < lowestCharge[0] and charge >= 0:
                lowestCharge = (charge, pH)
            pH += 0.01  # assigns pH to two decimals
        return (lowestCharge[1])

    def aaComposition(self):
        '''
        Creates an amino acids dictionary, and calculates amino acid composition and return it.
        '''
        compositionDict = {}  # creates empty dictionary for amino acid composition
        for amino in self.aa2mw.keys():  # access keys from dictionary
            # stores amino acids and composition values to dictionary
            compositionDict[amino] = self.proteinString.count(amino)
        return compositionDict

    def _charge_(self, pH):
        '''
        Calculates net charge in of protein from its sequence using the formula provided
        in the assigment description.
        '''
        # makes total charge equal to zero
        totalCharge = 0.00
        for amino in self.proteinString:
            # if amino acid is in the positive charge dictionary, it will be executed
            if amino in self.aa2chargePos.keys():
                powerA = 10 ** (self.aa2chargePos[amino])
                powerB = 10 ** (self.aa2chargePos[amino]) + 10 ** pH
                totalCharge += powerA / powerB
            # if amino acid is in the negative charge dictionary, it will be executed
            elif amino in self.aa2chargeNeg.keys():
                powerA = 10 ** pH
                powerB = 10 ** (self.aa2chargeNeg[amino]) + 10 ** pH
                totalCharge -= powerA / powerB
        # Calculates the Amino (N), and Carboxyl (C) terminus in protein sequence
        totalCharge += (10 ** self.aaNterm) / (
                    10 ** self.aaNterm + 10 ** pH)  # assigns value to total charge adding amino terminus
        totalCharge -= (10 ** pH) / (
                    10 ** self.aaCterm + 10 ** pH)  # assigns value to total charge by subtracting the carboxyl terminus
        pass

        return totalCharge

    def molarExtinction(self):
        '''
        Calcultes protein's molar extinction coefficient from amino acids compositions.
        Method adds the molar extinction coefficients of Y, W, and C at 280 nm absorbance found in
        protein sequence, and returns a float.
        '''
        # calculates molar extinction coefficient of tyrosine
        molarY = self.proteinString.count("Y") * self.aa2abs280["Y"]
        # calculates molar extinction coefficient of tryptophan
        molarW = self.proteinString.count("W") * self.aa2abs280["W"]
        # calculates molar extintion coefficient of cysteine
        molarC = self.proteinString.count("C") * self.aa2abs280["C"]
        return float(molarY + molarW + molarC)

    def massExtinction(self):
        '''
        Calculates mass extinction coefficient from the molar extinction coefficient
        '''
        molecularWeight = self.molecularWeight()
        return self.molarExtinction() / molecularWeight if molecularWeight else 0.0

    def molecularWeight(self):
        '''
        Calculates molecular weight of protein from its amino acids sequence
        '''
        molecularWeight = 0.0  # defines the weight of water
        # calculates the number of amino acids in protein sequence
        aaLength = sum(len(string) for string in self.proteinString)
        thisCount = self.mwH2O * (aaLength - 1)  # calcluates molecular weight
        for amino in self.proteinString:
            molecularWeight += self.aa2mw[amino]
        return (molecularWeight - thisCount)

#####################################################################################################
# Class FastAreader
#####################################################################################################

'''
In this class, we define objects to read FastaA files inside the sequenceAnalysis program.
'''

import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence


        