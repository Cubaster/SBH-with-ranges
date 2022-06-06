from random import randint, shuffle


class Generator:
    # Class generates n-length sequence of nitrogenous bases, divides sequence in k-length oligonucleotides
    # and shuffles them to create spectrum

    # ___ATTRIBUTES___
    # self.dnaLength - DNA sequence length
    # self.oligoLength - oligonucleotide length
    # self.oligoDict - dictionary holds number of oligonucleotides occurrence
    # self.starter - first oligonucleotide in sequence
    # self.__alphabet - holds nitrogenous bases: A - adenine, C - cytosine, G - guanine, T- thymine
    # self.sequence - DNA generated sequence

    # ___METHODS___
    # generateSequence(self) - generate DNA sequence
    # __getOligos(self) - create oligonucleotides from sequence
    # getSpectrum(self) - returns spectrum (shuffled oligonucleotides)

    def __init__(self, n: int, k: int):
        # ___PARAMS___
        # n - DNA sequence length
        # k - oligonucleotide length

        self.dnaLength = n
        self.oligoLength = k
        self.oligoDict = {}
        self.starter = None
        self.__alphabet = ['A', 'C', 'G', 'T']
        self.sequence = ''

    def generateSequence(self):
        # generate sequence by adding randomly alphabet elements

        for i in range(self.dnaLength):
            self.sequence += self.__alphabet[randint(0, 3)]
        self.__getOligos()

    def __getOligos(self):
        # create oligonucleotides from sequence
        # saves first oligonucleotide as starter
        # adds oligonucleotides to dictionary and counts oligos occurrence

        # ___Variables___
        # oligo - oligonucleotide is k-length substring of DNA sequence

        for i in range(self.dnaLength - self.oligoLength + 1):
            oligo = self.sequence[i:self.oligoLength + i]
            if i == 0:
                self.starter = oligo
            if oligo not in self.oligoDict:
                self.oligoDict[oligo] = 1
            else:
                self.oligoDict[oligo] += 1

    def getSpectrum(self):
        # returns shuffled oligonucleotides

        spectrum = list(self.oligoDict.keys())
        shuffle(spectrum)
        return spectrum


# TODO remove __main__ statement after class validation and verification


if __name__ == "__main__":
    generator = Generator(3, 4)
    generator.generateSequence()
    generator.getSpectrum()

