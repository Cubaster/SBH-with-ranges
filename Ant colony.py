from generator import Generator
from ant import Ant
from utilities import generatePheromonesMatrix, generateWeightsMatrix, createTranslation, mergeSolution
from copy import deepcopy


class AntColony:
    def __init__(self, ant_count: int = 50, alpha: float = 0.125, evaporation_coefficient: float = 0.4,
                 iterations: int = 80, sequence_length: int = 500, oligo_size: int = 9, percent: float = 0.05):
        self.percent = percent
        self.oligo_size = oligo_size
        self.sequence_length = sequence_length
        self.alpha = alpha
        self.ant_count = ant_count
        self.iterations = iterations
        self.evaporation_coefficient = evaporation_coefficient
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
        self.generator = Generator(self.sequence_length, self.oligo_size, self.percent)
        self.generator.generateSequence()
        self.starter = self.generator.starter
        self.initial_solution = self.generator.getSpectrum()
        self.ranges = self.generator.oligoDict

    def _initArea(self):
        self.weights = generateWeightsMatrix(self.sequence_length, self.initial_solution)
        self.translation = createTranslation(self.initial_solution)
        self.pheromones_map = generatePheromonesMatrix(self.sequence_length)

    def _init_ants(self):

        ant_list = []
        if self.first_attempt:
            self.first_attempt = False
            for ant in range(self.ant_count):
                new_ant = deepcopy(
                    Ant(self.starter, self.initial_solution, self.weights, self.translation, self.pheromones_map,
                        self.ranges, self.sequence_length, self.alpha, True))
                ant_list.append(new_ant)
            return ant_list

        for ant in range(self.ant_count):
            new_ant = deepcopy(
                Ant(self.starter, self.initial_solution, self.weights, self.translation, self.pheromones_map,
                    self.ranges, self.sequence_length, self.alpha))
            ant_list.append(new_ant)
        return ant_list

    def _updatePheromonesMap(self, pheromones_maps: list):
        for pheromones_map in pheromones_maps:
            for i in range(self.sequence_length):
                for j in range(self.sequence_length):
                    self.pheromones_map[i][j] += pheromones_map[i][j]

        for i in range(self.sequence_length):
            for j in range(self.sequence_length):
                self.pheromones_map[i][j] *= (1 - self.evaporation_coefficient)

    def _bestSolution(self, solution):
        new_solution = mergeSolution(solution, self.oligo_size)
        print(f"solution: {self.generator.sequence}")
        print(f"new solution: {new_solution}\n")

    def mainloop(self):
        self._initGenerator()
        self._initArea()
        self.ants = self._init_ants()
        for iteration in range(self.iterations):
            print(f"iteration {iteration}")
            pheromones_maps = []
            for ant in self.ants:
                print(ant.ranges)
                solution, pheromones_map = ant.run()
                pheromones_maps.append(pheromones_map)
                self._bestSolution(solution)
            self._updatePheromonesMap(pheromones_maps)
            self.ants = self._init_ants()


if __name__ == "__main__":
    ant_colony = AntColony(ant_count=10, iterations=20, sequence_length=200, oligo_size=5)
    ant_colony.mainloop()
