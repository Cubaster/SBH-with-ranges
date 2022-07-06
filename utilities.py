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


def generatePheromonesMatrix(size: int):
    """
    Create and initialize pheromones map to 0's
    :param size: size of sequence
    :return: pheromones map
    """
    pheromones_map = [[0.1 for _ in range(size)] for _ in range(size)]
    return pheromones_map


def generateWeightsMatrix(size: int, spectrum: list):
    """
    generate weights matrix base on oligonucleotides distance

    :param spectrum: list of nucleotide names
    :param size: size of sequence
    :return: weights matrix
    """

    weightsMatrix = [[0 for _ in range(size)] for _ in range(size)]
    length = len(spectrum)
    for i in range(length):
        for j in range(length):
            weightsMatrix[i][j] = customDistance(spectrum[i], spectrum[j])
    return weightsMatrix


def mergeSolution(solution: list, oligonucleotide_size: int):
    """
    Concatenate oligonucleotides in one sequence

    :param solution: solution created after ant journey
    :param oligonucleotide_size: oligonucleotide size
    :return: DNA sequence
    """
    result = solution[0]
    while len(solution) != 1:
        dist = customDistance(solution[0], solution[1])
        tmp = solution[1][oligonucleotide_size - dist: oligonucleotide_size]
        result += tmp
        solution.remove(solution[0])
    return result

