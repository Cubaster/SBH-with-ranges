from generator import Generator
from ant import Ant
from utilities import generatePheromonesMatrix, generateWeightsMatrix, createTranslation, mergeSolution, \
    levenshteinDistance, oligonucleotideComparison
from copy import deepcopy


class AntColony:
    def __init__(self, ant_count: int = 50, alpha: int = 2, evaporation_coefficient: float = 0.1,
                 iterations: int = 30, sequence_length: int = 700, oligo_size: int = 8, percent: float = 0.01, positive=False):
        """
        Initialize ant colony.

        :param ant_count: number of ants in colony
        :param alpha: value of trace left by ant (pheromones map)
        :param evaporation_coefficient: coefficient used in formula (1 - coeff) define trace evaporation
        :param iterations: number of journeys made by each ant
        :param sequence_length: DNA sequence length
        :param oligo_size: oligonucleotide size
        :param percent: define width of ranges in which oligonucleotide occurs
        """
        self.positive = positive
        self.percent = percent
        self.oligo_size = oligo_size
        self.sequence_length = sequence_length
        self.alpha = alpha
        self.ant_count = ant_count
        self.iterations = iterations
        self.evaporation_coefficient = evaporation_coefficient
        self.best_result = 10000000000000
        self.best_solution = None
        self.ants = None
        self.first_attempt = True

        self.generator = None
        self.starter = None
        self.initial_solution = None
        self.ranges = None

        self.weights = None
        self.translation = None
        self.pheromones_map = None

    def _initGenerator(self):
        """
        Initialize generator and sets information about sequence: sequence, starter, initial_solution(shuffled oligos), and ranges in which oligos occurs. .
        """
        self.generator = Generator(self.sequence_length, self.oligo_size, self.percent, self.positive)
        self.generator.generateSequence()
        self.starter = self.generator.starter
        self.initial_solution = self.generator.getSpectrum()
        self.ranges = self.generator.oligoDict

    def _initArea(self):
        """
        Initialize area (graph connections):
            weights (distances between oligos - cost),
            translation: map oligonucleotide to vertex number,
            pheromones_map: map of pheromone's trace (probabilities of choosing edge)
        """
        self.weights = generateWeightsMatrix(self.sequence_length, self.initial_solution)
        self.translation = createTranslation(self.initial_solution)
        self.pheromones_map = generatePheromonesMatrix(self.sequence_length, len(self.initial_solution))

    def _createAnt(self, firstAttempt: bool):
        """
        Create unique instance of ant
        :param firstAttempt: ant already started journey
        :return: ant instance11
        """
        starer = self.starter
        initial_solution = deepcopy(self.initial_solution)
        weights = deepcopy(self.weights)
        translation = self.translation
        pheromones_map = deepcopy(self.pheromones_map)
        ranges = deepcopy(self.ranges)
        sequence_length = self.sequence_length
        alpha = self.alpha
        oligo_size = self.oligo_size
        ant = Ant(starer, initial_solution, weights, translation, pheromones_map, ranges, sequence_length, alpha, oligo_size,
                  firstAttempt)
        return ant

    def _init_ants(self):
        """
        Initialize ants which will be set on journey
        :return: list of ants
        """
        ant_list = []
        if self.first_attempt:
            self.first_attempt = False
            for ant in range(self.ant_count):
                new_ant = self._createAnt(True)
                ant_list.append(new_ant)
            return ant_list

        for ant in range(self.ant_count):
            new_ant = self._createAnt(False)
            ant_list.append(new_ant)
        return ant_list

    def _updatePheromonesMap(self, pheromones_maps: list):
        """
        Update Pheromones rate (probability of choosing edge) after all ants return.

        :param pheromones_maps: list of pheromones maps for each ant
        """
        divide = 0
        # find norm (biggest element for this journey)
        for pheromones_map in pheromones_maps:
            for i in range(self.sequence_length):
                for j in range(self.sequence_length):
                    if pheromones_map[i][j] > divide:
                        divide = pheromones_map[i][j]

        # trace evaporation
        for i in range(self.sequence_length):
            for j in range(self.sequence_length):
                if self.pheromones_map[i][j] > 100:
                    self.pheromones_map[i][j] = 100
                else:
                    self.pheromones_map[i][j] += (pheromones_map[i][j] / divide) * (1 - self.evaporation_coefficient)

    def _bestSolution(self, solution, spots: list):
        """
        Finds best oligonucleotides

        :param solution: last sequence fount by ants
        """
        new_solution = mergeSolution(solution, self.oligo_size)
        result = levenshteinDistance(new_solution, self.generator.sequence)
        if result < self.best_result:
            self.best_result = result
            self.best_solution = new_solution

    def mainloop(self):
        """
        Initialize all actions in colony

        :return: best solution found by ants
        """
        self._initGenerator()
        self._initArea()
        self.ants = self._init_ants()
        pheromones_maps = []
        for iteration in range(self.iterations):
            print(f"iteration {iteration}")

            for ant in self.ants:
                solution, pheromones_map = ant.run()
                pheromones_maps.append(pheromones_map)
                self._bestSolution(solution, ant.cover_spots)
                del ant
            self._updatePheromonesMap(pheromones_maps)
            self.ants = self._init_ants()
        oligonucleotideComparison(self.generator.sequence, self.best_solution, self.oligo_size)
        return self.best_solution


if __name__ == "__main__":
    ant_count = 50
    iterations = 30
    sequence_length = 300
    oligo_size = 8
    alpha = 10
    evaporation_coeff = 0.1
    percent = 0.01
    positive = [True, False]

    ant_colony = AntColony(percent=percent, positive=False)
    best = ant_colony.mainloop()
    print(best)
    print(ant_colony.generator.sequence)
    print(ant_colony.best_result / sequence_length * 100)
