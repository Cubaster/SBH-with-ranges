from threading import Thread
from random import choice, choices


# callback_function triggers when length of path riches length od DNA sequence
# (Thread)
class Ant:
    def __init__(self, starting_point: str, nodes: list, weights: list, translations: dict, pheromones_map: list,
                 ranges_dict: dict, sequenceLength: int, alpha, first_attempt=False):
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
        # Thread.__init__(self)

        self.sequenceLength = sequenceLength
        self.starting_point = starting_point
        self.current_location = starting_point
        self.nodes = nodes
        self.weights = weights
        self.translations = translations
        self.pheromones_map = pheromones_map
        self.ranges = ranges_dict
        self.first_attempt = first_attempt
        self.route = []
        self.sequence_cover = 0

    def run(self):
        """
        ant perform actions as long as sequence_cover condition is not met
        """

        while self.sequence_cover <= self.sequenceLength:
            moveTo = self._choosePath()
            self._move(moveTo)

    def _verifyPath(self, node):
        """
        verify if node can be chosen - self.sequence between ranges

        :param node: oligonucleotide name e.g. "ATCGAT"
        :return: flag indicating if node can be chosen
        """
        for counter in range(0, len(self.ranges[node]), 2):
            if self.ranges[node][counter] <= self.sequence_cover <= self.ranges[node][counter + 1]:
                self.ranges[node].remove(self.ranges[node][counter + 1])
                self.ranges[node].remove(self.ranges[node][counter])
                return True, node
        return False, node

    @property
    def _choosePath(self):
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
        if self.first_attempt:
            while not go:
                go, node = self._verifyPath(choice(self.nodes))
            return node

        pheromonesForCurrentNode = self.pheromones_map[self.translations[self.current_location]]
        while not go:
            [oligonucleotide] = choices(self.nodes, pheromonesForCurrentNode)
            if oligonucleotide not in alreadyChecked:
                go, node = self._verifyPath(oligonucleotide)
                alreadyChecked.append(node)
        return node

    def _updatePath(self, moveTo):
        """
        add chosen node to route and remove oligonucleotide from nodes if there is no more occurrence

        :param moveTo: chosen node
        """
        self.route.append(moveTo)
        if not self.ranges[moveTo]:
            self.nodes.remove(moveTo)

    def _updateSequence(self, moveTo):
        self.sequence_cover += self.weights[self.translations[self.current_location]][self.translations[moveTo]]

    def _move(self, moveTo):
        self._updatePath(moveTo)
        self._updateSequence(moveTo)
        self.current_location = moveTo


