from asyncio.windows_events import NULL
from msilib import sequence
from unittest import result
from generator import Generator
from ant import Ant
from utilities import generatePheromonesMatrix, generateWeightsMatrix, createTranslation, mergeSolution, levenshteinDistance
from copy import deepcopy




class AntColony:
    def __init__(self, ant_count: int = 50, alpha: int = 12, evaporation_coefficient: float = 0.4,
                 iterations: int = 80, sequence_length: int = 500, oligo_size: int = 9, percent: float = 0.05):
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
        self.first_attempt = False

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

    def _createAnt(self, firstAttempt: bool):
        starer = self.starter
        initial_solution = deepcopy(self.initial_solution)
        weights = deepcopy(self.weights)
        translation = self.translation
        pheromones_map = deepcopy(self.pheromones_map)
        ranges = deepcopy(self.ranges)
        sequence_length = self.sequence_length
        alpha = self.alpha
        ant = Ant(starer, initial_solution, weights, translation, pheromones_map, ranges, sequence_length, alpha, firstAttempt)
        return ant

    def _init_ants(self):

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
        for pheromones_map in pheromones_maps:
            for i in range(self.sequence_length):
                for j in range(self.sequence_length):
                    if self.pheromones_map[i][j] > 100:
                        self.pheromones_map[i][j] = 100
                    else:
                        self.pheromones_map[i][j] += pheromones_map[i][j]/self.ant_count

        for i in range(self.sequence_length):
            for j in range(self.sequence_length):
                self.pheromones_map[i][j] *= (1 - self.evaporation_coefficient)

    def _bestSolution(self, solution):
        new_solution = mergeSolution(solution, self.oligo_size)
        #print(f"solution: {self.generator.sequence}")
        #print(f"new solution: {new_solution}\n")
        result = levenshteinDistance(new_solution, self.generator.sequence)
        #print(result)
        if result < self.best_result:
             self.best_result = result
             self.best_solution = new_solution
        
    def mainloop(self):
        self._initGenerator()
        self._initArea()
        self.ants = self._init_ants()
        pheromones_maps = []
        for iteration in range(self.iterations):
            #print(f"iteration {iteration}")
            
            #print(pheromones_maps)
            for ant in self.ants:
                solution, pheromones_map = ant.run()
                pheromones_maps.append(pheromones_map)
                self._bestSolution(solution)
            self._updatePheromonesMap(pheromones_maps)
            self.ants = self._init_ants()
        return self.best_solution, self.initial_solution


if __name__ == "__main__":
    ant_count = [1,5,10,20]
    iterations = [1,3,5,7]
    sequence_length = 30; oligo_size = 4
    alpha = 100; evaporation_coeff = 0.2; percent = 0.1
    for count in ant_count:
        for iter in iterations:
            print(f"Ants: {count} iterations: {iter}")
            ant_colony = AntColony(ant_count = count, alpha=alpha, evaporation_coefficient=evaporation_coeff, iterations = iter, sequence_length = sequence_length, oligo_size = oligo_size)
            best, initial = ant_colony.mainloop()
            print(ant_colony.best_result /sequence_length*100)
            ant_colony = NULL
            #print(initial)
            #print(best)
            #print(ant_colony.pheromones_map)
    

