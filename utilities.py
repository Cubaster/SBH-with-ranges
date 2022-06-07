
def customDistance(oligo1: str, oligo2: str):
    """
    return index no. at which oligo1 substring is identical with oligo2 substring - further interpreted as weight
    strings has to be equal in length

    :param oligo1: first oligonucleotide
    :param oligo2: second oligonucleotide
    :return: weight of similarity - the bigger weight the smaller similarity
    """
    #

    rowSize = len(oligo1)
    columnSize = len(oligo2)
    if rowSize != columnSize:
        return None

    for i in range(1, rowSize):
        if oligo1[i:rowSize] == oligo2[0: rowSize - i]:
            return i
    return rowSize


def createTranslation(oligonucleotides: list):
    """
    translate oligonucleotide names to numeric values corresponding to matrix weights

    :param oligonucleotides: spectrum produced by generator
    :return: dictionary holding fields number in matrix of weights
    """

    translation = {}
    counter = 0
    for oligonucleotide in oligonucleotides:
        translation[oligonucleotide] = counter
        counter += 1

    return translation




