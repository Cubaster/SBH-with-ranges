
def customDistance(oligo1: str, oligo2: str):
    # return index no. at which oligo1 whole string\string end is identical with oligo2 whole string\substring beginning
    # strings has to be equal in length

    rowSize = len(oligo1)
    columnSize = len(oligo2)
    if rowSize != columnSize:
        return None

    for i in range(rowSize):
        if oligo1[i:rowSize] == oligo2[0: rowSize - i]:
            return i
    return rowSize


