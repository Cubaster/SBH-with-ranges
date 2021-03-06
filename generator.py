from random import randint, shuffle, choices


class Generator:
    """
    Class generates n-length sequence of nitrogenous bases, divides sequence in k-length oligonucleotides
    and shuffles them to create spectrum

    ___ATTRIBUTES___
    self.dnaLength - DNA sequence length
    self.oligoLength - oligonucleotide length
    self.oligoDict - dictionary holds number of oligonucleotides occurrence
    self.starter - first oligonucleotide in sequence
    self.percent - percent of sequence which extends range of oligonucleotide position
    self._alphabet - holds nitrogenous bases: A - adenine, C - cytosine, G - guanine, T- thymine
    self.sequence - DNA generated sequence
    self.positive - generating positive mistakes

     ___METHODS___
    generateSequence(self) - generate DNA sequence
    _getOligonucleotides(self) - create oligonucleotides from sequence
    getSpectrum(self) - returns spectrum (shuffled oligonucleotides)
    _getRange(self, position) - return range in which oligonucleotide is stored
    _addToDict - add oligonucleotides(ranges) to dictionary
    _generatePositive(self) - generate positive instances (not existing oligonucleotides)
    """

    def __init__(self, n: int, k: int, percent: float, positive=False):
        """
        generator initialization

        :param n: DNA sequence length
        :param k: oligonucleotide length
        :param percent: percent of extra range to basic position
        :param positive: flag decides whether generate positive mistakes
        """

        self.dnaLength = n
        self.oligoLength = k
        self.oligoDict = {}
        self.starter = None
        self.percent = percent
        self._alphabet = ['A', 'C', 'G', 'T']
        self.sequence = ''
        self.positive = positive

    def generateSequence(self):
        """
        generate sequence by adding random alphabet elements
        """

        for i in range(self.dnaLength):
            self.sequence += self._alphabet[randint(0, 3)]
        self._getOligonucleotides()

    def _getOligonucleotides(self):
        """
        create oligonucleotides from sequence
        saves first oligonucleotide as starter

        :parameter oligo: oligonucleotide is k-length substring of DNA sequence
        """

        for i in range(self.dnaLength - self.oligoLength + 1):
            oligo = self.sequence[i:self.oligoLength + i]
            if i == 0:
                self.starter = oligo
            self._addToDict(oligo, i)

        if self.positive:
            self._generatePositive()

    def _addToDict(self, oligo, position):
        """
         add oligonucleotides to dictionary and append its lower and upper range to flat list
        e.g. [2, 4, 67, 81] lowers: 2, 67 uppers: 4, 81

        :param oligo: oligonucleotide
        :param position: oligonucleotide exact position in DNA sequence
        """
        if oligo not in self.oligoDict:
            self.oligoDict[oligo] = []
            low, high = self._getRange(position)
            self.oligoDict[oligo].append(low)
            self.oligoDict[oligo].append(high)
        else:
            low, high = self._getRange(position)
            self.oligoDict[oligo].append(low)
            self.oligoDict[oligo].append(high)

    def getSpectrum(self):
        """
        shuffle oligonucleotides dictionary keys

        :return: set of dictionary keys (oligonucleotides)
        """
        #

        spectrum = list(self.oligoDict.keys())
        shuffle(spectrum)
        return spectrum

    def _getRange(self, position):
        """
        generate range in which oligonucleotide should be placed during reconstruction

        :param position: exact position of nucleotide in DNA sequence
        :return: range in which oligonucleotide should be placed during reconstruction
        """
        lowerRange = position - randint(0, int(self.dnaLength * self.percent))
        if lowerRange < 0:
            lowerRange = 0
        upperRange = position + randint(0, int(self.dnaLength * self.percent))
        if upperRange > self.dnaLength:
            upperRange = self.dnaLength
        oligoRange = [lowerRange, upperRange]
        return oligoRange

    def _generatePositive(self):
        """
        create and add to dict artificial instances on nucleotides
        active when :param: positive is set to True
        """
        iterations = int(self.dnaLength * 0.1)
        for iteration in range(iterations):
            oligonucleotide = ''
            for base in range(self.oligoLength):
                oligonucleotide += self._alphabet[randint(0, 3)]
            self._addToDict(oligonucleotide, randint(0, self.dnaLength))

    # TODO remove __main__ statement after class validation and verification

if __name__ == "__main__":
    generator = Generator(70, 3, 0.1, True)
    generator.generateSequence()
    # print(generator.sequence)
    #print(generator.getSpectrum())
