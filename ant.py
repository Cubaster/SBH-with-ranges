from copy import deepcopy
from random import choice, choices


class Ant:
    def __init__(self, starting_point: str, nodes: list, weights: list, translations: dict, pheromones_map: list,
                 ranges_dict: dict, sequenceLength: int, alpha: float, oligo_size: int, first_attempt=False):
        """
        Initialize ant graph traverse

        :param starting_point: ant starting point e.g "ATCG"
        :param sequenceLength: length of DNA sequence
        :param nodes: nodes which can be visited - during journey its number decrease (visited nodes) e.g. ["ACT", "AGT" ..]
        :param weights: matrix of weights - value of difference between oligonucleotides e.g.[[1, 3, 4][1, 2, 4][7, 10, 1]]
        :param translations: dictionary binding nodes names(oligonucleotides names) with weights e.g. {"ACT": 0, "AGT": 1}
        :param pheromones_map: pheromone trace to change choice probability e.g. [[0,0,0.1][0,0.2,0][1,0.5,0.2]]
        :param ranges_dict: holds ranges defining when some node should be visited e.g{"ACT": [0, 3, 7, 13], "AGT": [4, 7]}
        :param alpha: pheromone rate
        :param oligo_size: oligonucleotide size
        :param first_attempt: when true random choices are enabled

        :parameter route: store visited nodes
        :parameter stop_walk: flag indicating if traversal is finished
        :parameter sequence_cover: sum of weights during traversal if reaches limit (value of DNA length) ant stops
        """

        self.oligo_size = oligo_size
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
        self.cover_spots = []

    def run(self):
        """
        Ant perform actions as long as sequence_cover condition is not met: choose path and move to next vertex
        """
        self._move(self.starting_point)
        while self.sequence_cover < self.sequenceLength:
            moveTo = self._choosePath()
            self._move(moveTo)
        return self.route, self.pheromones_map

    def _verifyPath(self, node):
        """
        Verify if node can be chosen - current value within ranges (valid): remove ranges of used vertex

        :param node: oligonucleotide name e.g. "ATCGAT"
        :return: flag indicating if node can be chosen and node itself
        """

        for counter in range(0, len(self.ranges[node]), 2):
            if self.ranges[node][counter] <= self.sequence_cover <= self.ranges[node][counter + 1]:
                self.ranges[node].remove(self.ranges[node][counter + 1])
                self.ranges[node].remove(self.ranges[node][counter])
                return True, node
        return False, node

    def _getWeights(self, nodes: list):
        """
        Function get pheromone trace between current node and possible nodes.

        :param nodes: list on yet not used nodes
        :return: list of pheromones
        """
        base_list = self.pheromones_map[self.translations[self.current_location]]
        index_list = [self.translations[element] for element in nodes]
        final_list = []
        for index in index_list:
            final_list.append(base_list[index])
        return final_list

    def _choosePath(self):
        """
        During first attempt algorithm search for node to which ant can go at random

        Each other attempt consider k elements base on pheromones and choose best fitting.

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

        # first attempt choose at random
        if self.first_attempt:
            while not go:
                oligonucleotide = choice(self.nodes)
                if oligonucleotide not in alreadyChecked:
                    go, node = self._verifyPath(oligonucleotide)
                    alreadyChecked.append(node)
            return node

        pheromonesForCurrentNode = self._getWeights(nodes)
        while not go:
            # chose k oligonucleotides base on weights (random choice)
            oligonucleotides = choices(nodes, weights=pheromonesForCurrentNode, k=self.oligo_size)
            best_oligo = self.oligo_size + 1
            oligonucleotide = ''
            for oligos in oligonucleotides:
                if self.weights[self.translations[self.current_location]][self.translations[oligos]] < best_oligo:
                    best_oligo = self.weights[self.translations[self.current_location]][self.translations[oligos]]
                    oligonucleotide = oligos
            # reject bad fitting oligonucleotides
            if best_oligo < 5:
                if oligonucleotide not in alreadyChecked and oligonucleotide != '':
                    go, node = self._verifyPath(oligonucleotide)
                    alreadyChecked.append(node)
            # choose last tried if any do not meet condition
            if counter > 3 * self.sequenceLength:
                return oligonucleotide
            counter += 1
        return node

    def _updatePath(self, moveTo):
        """
        Add chosen node to route and remove oligonucleotide from nodes if there is no more occurrence

        :param moveTo: chosen node
        """
        self.route.append(moveTo)
        if not self.ranges[moveTo]:
            if moveTo in self.nodes:
                self.nodes.remove(moveTo)
            else:
                print("positive")

    def _updateSequence(self, moveTo):
        """
        Update sequence cover base on distance between nodes and mark path with pheromones
        :param moveTo: chosen oligonucleotide
        """
        self.sequence_cover += self.weights[self.translations[self.current_location]][self.translations[moveTo]]
        self.cover_spots.append(self.sequence_cover)
        self.pheromones_map[self.translations[self.current_location]][self.translations[moveTo]] += self.alpha

    def _move(self, moveTo):
        """
        Perform ant movement action

        :param moveTo: chosen oligonucleotide
        """
        self._updatePath(moveTo)
        self._updateSequence(moveTo)
        self.current_location = moveTo
