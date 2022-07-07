from copy import deepcopy
from email import generator
from threading import Thread
from random import choice, choices
from generator import Generator
from utilities import *
from time import sleep


# callback_function triggers when length of path riches length od DNA sequence
# (Thread)
class Ant:
    def __init__(self, starting_point: str, nodes: list, weights: list, translations: dict, pheromones_map: list,
                 ranges_dict: dict, sequenceLength: int, alpha: float, first_attempt=False):
        """
        initialize ant graph traverse

        :param starting_point: ant starting point e.g "ATCG"
        :param sequenceLength: length of DNA sequence
        :param nodes: nodes which can be visited - during journey its number decrease (visited nodes) e.g. ["ACT", "AGT" ..]
        :param weights: matrix of weights - value of difference between oligonucleotides e.g.[[1, 3, 4][1, 2, 4][7, 10, 1]]
        :param translations: dictionary binding nodes names(oligonucleotides names) with weights e.g. {"ACT": 0, "AGT": 1}
        :param pheromones_map: pheromone trace to change choice probability e.g. [[0,0,0.1][0,0.2,0][1,0.5,0.2]]
        :param ranges_dict: holds ranges defining when some node should be visited e.g{"ACT": [0, 3, 7, 13], "AGT": [4, 7]}
        :param alpha: pheromone rate
        :param first_attempt: when true random choices are enabled

        :parameter route: store visited nodes
        :parameter stop_walk: flag indicating if traversal is finished
        :parameter sequence_cover: sum of weights during traversal if reaches limit (value of DNA length) ant stops
        """
        #Thread.__init__(self)

        self.sequenceLength = sequenceLength
        self.starting_point = starting_point
        self.current_location = starting_point
        self.nodes = nodes
        self.weights = weights
        self.translations = translations
        self.pheromones_map = pheromones_map
        self.ranges = ranges_dict
        self.first_attempt = first_attempt
        self.alpha = alpha
        self.route = []
        self.sequence_cover = 0

    def run(self):
        """
        ant perform actions as long as sequence_cover condition is not met
        """
        #self._updatePath(self.starting_point)
        self._move(self.starting_point)
        while self.sequence_cover < self.sequenceLength:
            moveTo = self._choosePath()
            self._move(moveTo)

        return self.route, self.pheromones_map

    def _verifyPath(self, node):
        """
        verify if node can be chosen - self.sequence between ranges

        :param node: oligonucleotide name e.g. "ATCGAT"
        :return: flag indicating if node can be chosen
        """

        for counter in range(0, len(self.ranges[node]), 2):
            # print(f"sequence cover: {self.sequence_cover}")
            # print(f"range: {self.ranges[node][counter]} {self.ranges[node][counter+1]}")
            # print(self.ranges[node][counter] <= self.sequence_cover <= self.ranges[node][counter + 1])
            if self.ranges[node][counter] <= self.sequence_cover <= self.ranges[node][counter + 1]:
                self.ranges[node].remove(self.ranges[node][counter + 1])
                self.ranges[node].remove(self.ranges[node][counter])
                # print("tutaj")
                return True, node
        return False, node

    def _getWeights(self, nodes: list):
        base_list = self.pheromones_map[self.translations[self.current_location]]
        index_list = [self.translations[element] for element in nodes]
        final_list = []
        for index in index_list:
            final_list.append(base_list[index])
        return final_list

    def _choosePath(self):
        # print("boom!")
        """
        During first attempt algorithm search for node to which ant can go at random

        Each other attempt choose random path but choice is base on pheromone weights

        :var node: represents oligonucleotide name
        :var alreadyChecked: list of rejected nucleotides
        :var go: flag indicating finding next move
        :var pheromonesForCurrentNode: list of pheromone weights
        :var oligonucleotide: chosen oligonucleotide

        :return: node name
        """
        node = None
        alreadyChecked = []
        go = False
        counter = 0
        nodes = deepcopy(self.nodes)

        if self.first_attempt:
            # print("first")
            while not go:
                #sleep(1)
                oligonucleotide = choice(self.nodes)
                #print(self.nodes)
                # print(oligonucleotide)
                # print(self.ranges)
                # print(alreadyChecked)
                if oligonucleotide not in alreadyChecked:
                    go, node = self._verifyPath(oligonucleotide)
                    #print(go)
                    alreadyChecked.append(node)
                    
            return node
        #print(f"cover after: {self.sequence_cover}")
        pheromonesForCurrentNode = self._getWeights(nodes)
        # print(pheromonesForCurrentNode)
        while not go:
            #print("second")
            [oligonucleotide] = choices(nodes, weights=pheromonesForCurrentNode, k=1)
            # print(oligonucleotide)
            if oligonucleotide not in alreadyChecked:
                # print(oligonucleotide)
                # print(alreadyChecked)
                go, node = self._verifyPath(oligonucleotide)
                #print(go)
                alreadyChecked.append(node)
            if counter > self.sequenceLength:
                return oligonucleotide
            counter += 1
        return node

    def _updatePath(self, moveTo):
        """
        add chosen node to route and remove oligonucleotide from nodes if there is no more occurrence

        :param moveTo: chosen node
        """
        self.route.append(moveTo)
        if not self.ranges[moveTo]:
            if moveTo in self.nodes:
                self.nodes.remove(moveTo)
            else:
                print("positive")
        #print(self.ranges)

    def _updateSequence(self, moveTo):
        """
        Update sequence cover base on distance between nodes and mark path with pheromones
        :param moveTo: chosen oligonucleotide
        """
        self.sequence_cover += self.weights[self.translations[self.current_location]][self.translations[moveTo]]
        # print(self.sequence_cover)
        self.pheromones_map[self.translations[self.current_location]][self.translations[moveTo]] += self.alpha

    def _move(self, moveTo):
        """
        Perform general ant movement action

        :param moveTo: chosen oligonucleotide
        """
        self._updatePath(moveTo)
        self._updateSequence(moveTo)
        self.current_location = moveTo
        
# if __name__ == "__main__":
#     sequence_length = 70
#     oligonucleotide_lenght = 4
#     generator = Generator(sequence_length, oligonucleotide_lenght, 0.4)
#     generator.generateSequence()
#     spectrum = generator.getSpectrum()
#     #print(generator.sequence)
#     #print(spectrum)
#     weights = generateWeightsMatrix(70, spectrum)
#     translation = createTranslation(spectrum)
#     phermones = generatePheromonesMatrix(70)
#     ranges = generator.oligoDict
#     #print(generator.oligoDict)
#     # ant0 = Ant(generator.starter, newSpectrum(spectrum), weights, translation, phermones, newRanges(ranges), sequence_length, 0.25, True)
#     # ant1 = Ant(generator.starter, newSpectrum(spectrum), weights, translation, phermones, newRanges(ranges), sequence_length, 0.25)
#     # ant2 = Ant(generator.starter, newSpectrum(spectrum), weights, translation, phermones, newRanges(ranges), sequence_length, 0.25)
#     # ants = [ant0, ant1, ant2]
#     # for i in range(3):
#     #     print(f"Run {i}")
#     #     route, _ = ants[i].run()
#     #     print(ants[i].nodes)
#     #     print(route)
#     #     print(mergeSolution(route, oligonucleotide_lenght))
