#!/usr/bin/env python3

#####################################################################################################
# 
#   Student: Carlos Arevalo (caeareva)
#   Group: None
#
#   File: findORFs.py 
#   One fasta execution: python findORFs.py -lG -s "ATG" -mG 300 < tass2.fa > tass2ORFdata-ATG-100.txt
#   Multiple fastas execution: python findORFs.py -lG -s "ATG" -mG 0 < lab5test.fa > tass2ORFdata-ATG-100.txt
#   Pupose: find open reading frames in the complement and reverse complement of a fasta file.
#           Program was built to be executed in stdin and stdout.
#
#
#####################################################################################################

import sequenceAnalysis
import sys

class CommandLine():
    '''
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.
    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''

    def __init__(self, inOpts=None):
        '''
        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='findORFs.py - finds ORFs in a fasta file',
            epilog='Program epilog - please provide corrections and implementations for the program',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store', nargs='?', const=True, default=True,
                                 help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=range(0, 1000), default=100, action='store',
                                 help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action='append', nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

#####################################################################################################
# ORFfinder Class
#####################################################################################################

class ORFfinder():
    '''
    This class find open reading frames (ORFs) in an input fasta file. An ORF is defined 
    from the position of the start codon to the position of a stop codon. First, the class 
    finds ORFs in the complement of the double-stranded DNA sequence. Then, it finds the 
    ORFs in the reverse complement strand. Each strand has three possibe frames translations 
    both strands have a total of 6 possible ORFs: -3, -2, -1, 1, 2, 3. Method uses the 
    FastAreader class to read over the input fasta file, and then, it uses the ORFfinder class
    to find ORFs.
    '''
    complement = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'} # DNA complement dictionary
    orfsList = [[], [], []] # creates a list of list
    startPosition = [] # list stores the found start codon positions  
    stopPosition = [] # list stores the stop codons positions
    
    def __init__(self, seq):
        '''Initialize the program and create list for stop and start codons'''
        self.seq = seq
        self.inSeq = seq.replace(' ', '') # removes spaces in fasta sequence
        self.startCodon = ['ATG'] # start codon list 
        self.stopCodon = ['TAG', 'TAA', 'TGA'] # stop codons list

    def findORF(self):
        '''
        Find ORFs on complement strand and return a list. Define a codon as 3 nucleotides
        and find the positions of stop and start codons for each frame: 1, 2, 3.
        '''
        for frame in range(3):
            #self.stopPosition.clear()
            for position in range(frame, len(self.inSeq), 3):
                # defines a codon as 3 nucleotides 
                codon = self.inSeq[position: position+3]
                # if codon is a start codon save its position to list
                if codon in self.startCodon:
                    self.startPosition.append(position)
                    pass

                 # if codon is a stop codon save its positio to list
                if codon in self.stopCodon:
                    self.stopPosition.append(position)
                    # save frame and position if start position was not found and position equal to 1
                    if self.startPosition:
                        if not self.orfsList[frame] and len(self.stopPosition) == 1:
                            # save ORFs elements if found position to list
                            self.saveORF(0, position+3, position+3, frame)
                            self.startPosition.clear()
                        else: # if position not found, save ORFs elements
                            length = (position + 3) - self.startPosition[0]
                            self.saveORF(self.startPosition[0], (position+3), length,frame)
                            self.startPosition.clear()
                            pass

                    # save frame and position if start position was not found and position equal to 1
                    if not self.orfsList[frame] and len(self.stopPosition) == 1:
                        self.saveORF(0, position+3, position+3, frame)

                # save length and frame when position and start codon found
                if position == len(self.inSeq) - 4:
                    if self.startPosition:
                        length = (len(self.inSeq) - 1) - self.startPosition[0]
                        self.saveORF(self.startPosition[0], (len(self.inSeq) - 1), length, frame)
                        pass

        return self.orfsList

    def findReverseORF(self):
        '''
        Find ORFs on the reverse complement strand and return a list. Method reuses the findORF 
        function to find ORFs in the reverse complement strand: -1, -2, -3.
        '''
        self.orfsList = [[], [], []] # creates a list of list
        self.startPosition.clear() # clears start codon positions list 
        self.stopPosition.clear() # clears stop codons positions list
        pass
        
        # reverse complement sequence 
        reverseComp = self.reverseComplement() 
        # find ORFs in the reverse complement
        self.inSeq = reverseComp  
        return self.findORF()

    def saveORF(self, start, stop, length, frame):
        '''
        Save ORFs found in the complement and reverse complement strands in a list of list
        '''
        self.orfsList[frame].append((frame, start, stop, length))        

    def reverseComplement(self): 
        '''
        Define and return reverse complement of input DNA sequence
        '''
        return ''.join(self.complement.get(base, base) for base in self.inSeq[::-1])

def main(inCL=None):
    '''
    Find some genes. 
    Read in fasta file in stdin, calculates ORFs: frames, starts, stops, and lengths, 
    and returns a file. Function calculates the ORFs in complement strand first, and then,
    it calculates the ORFs in reverse complement strand and uses the stdout method to return
    the output file with ORFs.
    '''
    framesList = [] 
    if inCL is None:
        myCommandLine = CommandLine()
        if myCommandLine.args.longestGene:
            fastaFile = sequenceAnalysis.FastAreader()
            framesList.clear()
            # reads fasta file 
            for header, sequence in fastaFile.readFasta():
                # print header in stdout format
                sys.stdout.write(header + '\n')
                # read sequence in fasta and call class
                myFinder = ORFfinder(sequence)
                # find ORFs in complement strand 
                orfList = myFinder.findORF()
                # find ORFs in reverse complement strand
                reverseFrames = myFinder.findReverseORF()
                pass

                # acces ORFs in the complement strand
                for list in orfList:
                    # define elements complement in ORFs
                    for element in list:
                        frame = element[0] + 1 # find ORFs frame
                        start = element[1] + 1 # find ORFs start position
                        stop = element[2] # find ORFs stop position
                        # append all ORFs elements to list
                        framesList.append((frame, start, stop, element[3]))
                        pass

                # access ORFs in the reverse complement strand
                for list in reverseFrames:
                    # define elements in reverse ORFs
                    for element in list:
                        frame = element[0] + 1 # find ORFs' frame
                        start = len(sequence) - (element[2]) + 1 # find ORFs' start position
                        stop = len(sequence) - (element[1] + 1)  + 1 # find ORFs' stop position
                        # append all reverse ORFs' elements to list
                        framesList.append((-frame, start, stop, element[3]))
                        pass

                # sort ORFs and print file
                framesList.sort(key=lambda tup: (tup[3], tup[1]), reverse=True)
                for orf in framesList:
                    # returns or greater or equal to the minum gene length argument called
                    if orf[3] >= myCommandLine.args.minGene:
                        sys.stdout.write('{:+d} {:>5d}..{:>5d} {:>5d}\n'.format(orf[0], orf[1], orf[2], orf[3]))
                    
    else:
        myCommandLine = CommandLine(inCL)

if __name__ == "__main__":
    main()



