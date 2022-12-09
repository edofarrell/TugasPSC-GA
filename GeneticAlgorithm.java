
import java.util.Arrays;
import java.util.Comparator;

/**
 *
 * @author 
 *  Kelompok 15: 
 *      Keannen Renaldo Halim   - 6182001007
 *      Neil Christopher        - 6182001010 
 *      Edo Farrell Haryanto    - 6182001025
 */

public class GeneticAlgorithm {

    private int n;                  //ukuran papan (nxn)
    private int[][] board;          //array 2d unutk menyimpan papan permainan
    
    private int numOfGeneration;    //banyak generasi
    private int populationSize;     //besar populasi
    private double mutationRate;    //probabilitas terjadi mutasi
    private double crossoverRate;   //probabilitas crossover berhasil
    private int elitismCount;       //jumlah individu yang akan dipilih secara elitism

    //constructor
    public GeneticAlgorithm(int[][] board, int numOfGeneration, int populationSize, double mutationRate, double crossoverRate, int elitismCount) {
        this.n = board.length;                  //panjang papan permainan
        this.board = board;                     //papan permainan
        this.numOfGeneration = numOfGeneration; //banyak generasi
        this.populationSize = populationSize;   //besar populasi
        this.mutationRate = mutationRate;       //probabilitas terjadi mutasi
        this.crossoverRate = crossoverRate;     //probabilitas crossover berhasil
        this.elitismCount = elitismCount;       //jumlah individu yang akan dipilih secara elitism
    }

    //method untuk inisialisasi populasi
    public Population initPopulation(int chromosomeLength) {
        //inisialisasi populasi sesuai ukuran populasi dan panjang kromosom yang ditentukan
        Population population = new Population(this.populationSize, chromosomeLength);
        return population;
    }

    /**
     * Calculate fitness for an individual.
     *
     * In this case, the fitness score is very simple: it's the number of ones
     * in the chromosome. Don't forget that this method, and this whole
     * GeneticAlgorithm class, is meant to solve the problem in the "AllOnesGA"
     * class and example. For different problems, you'll need to create a
     * different version of this method to appropriately calculate the fitness
     * of an individual.
     *
     * @param individual the individual to evaluate
     * @return double The fitness value for individual
     */
    public double calcFitness(Individual individual) {

//		// Track number of correct genes
//		int correctGenes = 0;
//
//		// Loop over individual's genes
//		for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
//			// Add one fitness point for each "1" found
//			if (individual.getGene(geneIndex) == 1) {
//				correctGenes += 1;
//			}
//		}
//
//		// Calculate fitness
//		double fitness = (double) correctGenes / individual.getChromosomeLength();
//
//		// Store fitness
//		individual.setFitness(fitness);
//
//		return fitness;
        double fitness = 0;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (board[i][j] >= 0) {
                    int correctGenes = countNeighbour(i, j, board[i][j], individual);
                    if (correctGenes != -1) {
                        fitness += correctGenes;
                    } else {
                        fitness = 0;
                        individual.setFitness(fitness);
                        return fitness;
                    }
                }
            }
        }

        individual.setFitness(fitness);
        return fitness;
    }

    private int countNeighbour(int i, int j, int target, Individual individual) {
        int[][] move = {
            {-1, -1},
            {-1, 0},
            {-1, 1},
            {0, -1},
            {0, 0},
            {0, 1},
            {1, -1},
            {1, 0},
            {1, 1}
        };

        int count = 0;
        for (int k = 0; k < move.length; k++) {
            int newI = i + move[k][0];
            int newJ = j + move[k][1];
            if (validate(newI, newJ)) {
                int offset = n * newI + newJ;
                if (individual.getGene(offset) == 1) {
                    count++;
                }
                if (count > target) {
                    count = -1;
                    break;
                }
            }
        }

        return count;
    }

    private boolean validate(int i, int j) {
        return i < n && i >= 0 && j < n && j >= 0;
    }

    /**
     * Evaluate the whole population
     *
     * Essentially, loop over the individuals in the population, calculate the
     * fitness for each, and then calculate the entire population's fitness. The
     * population's fitness may or may not be important, but what is important
     * here is making sure that each individual gets evaluated.
     *
     * @param population the population to evaluate
     */
    public void evalPopulation(Population population) {
        double populationFitness = 0;

        // Loop over population evaluating individuals and suming population
        // fitness
        for (Individual individual : population.getIndividuals()) {
            populationFitness += calcFitness(individual);
            // populationFitness = Math.max(populationFitness, calcFitness(individual));
        }

        population.setPopulationFitness(populationFitness);

//        System.out.println("Next Generation:");
//        for (int i = 0; i < this.populationSize; i++) {
//            System.out.println(population.getIndividual(i).toString() + " " + population.getIndividual(i).getFitness());
//        }
//        System.out.println("");
    }

    /**
     * Check if population has met termination condition
     *
     * For this simple problem, we know what a perfect solution looks like, so
     * we can simply stop evolving once we've reached a fitness of one.
     *
     * @param population
     * @return boolean True if termination condition met, otherwise, false
     */
//    public boolean isTerminationConditionMet(Population population) {
//        for (Individual individual : population.getIndividuals()) {
//            if (individual.getFitness() == 1) {
//                return true;
//            }
//        }
//
//        return false;
//    }
    public boolean isTerminationConditionMet(int generation) {
        if (generation < this.numOfGeneration) {
            return false;
        } else {
            return true;
        }
    }

    /**
     * Select parent for crossover
     *
     * @param population The population to select parent from
     * @return The individual selected as a parent
     */
    public Individual selectParentRoulette(Population population) {
        // Get individuals
        Individual individuals[] = population.getIndividuals();

        // Spin roulette wheel
        double populationFitness = population.getPopulationFitness();
        double rouletteWheelPosition = Math.random() * populationFitness;

        // Find parent
        double spinWheel = 0;
        for (Individual individual : individuals) {
            spinWheel += individual.getFitness();
            if (spinWheel >= rouletteWheelPosition) {
                return individual;
            }
        }
        return individuals[population.size() - 1];
    }

    public Individual selectParentRank(Population population) {
        // Get individuals
        Individual individuals[] = population.getIndividuals();
        Arrays.sort(individuals, new Comparator<Individual>() {
            @Override
            public int compare(Individual o1, Individual o2) {
                if (o1.getFitness() > o2.getFitness()) {
                    return -1;
                } else if (o1.getFitness() < o2.getFitness()) {
                    return 1;
                }
                return 0;
            }
        });

        double totalRank = (individuals.length / 2.0) * (1 + individuals.length);
        double rank[] = new double[individuals.length];
        for (int i = individuals.length - 1; i >= 0; i--) {
            rank[i] = 1.0 * (individuals.length - i) / totalRank;
        }

        // Spin roulette wheel
        double rouletteWheelPosition = Math.random();

        // Find parent
        double spinWheel = 0;
        for (int i = 0; i < rank.length; i++) {
            spinWheel += rank[i];
            if (spinWheel >= rouletteWheelPosition) {
                return individuals[i];
            }
        }

        return individuals[population.size() - 1];
    }

    public Individual selectParentTournament(Population population) {
        // Get individuals
        Individual individuals[] = population.getIndividuals();
        Arrays.sort(individuals, new Comparator<Individual>() {
            @Override
            public int compare(Individual o1, Individual o2) {
                if (o1.getFitness() > o2.getFitness()) {
                    return -1;
                } else if (o1.getFitness() < o2.getFitness()) {
                    return 1;
                }
                return 0;
            }
        });

        int memberPerGroup = 5;
        double newPopulationFitness = 0;

        Individual newIndividuals[] = new Individual[individuals.length / memberPerGroup];
        for (int i = 0; i < newIndividuals.length; i++) {
            newIndividuals[i] = individuals[i * memberPerGroup];
            newPopulationFitness += newIndividuals[i].getFitness();
        }

        // Spin roulette wheel
        double rouletteWheelPosition = Math.random() * newPopulationFitness;

        // Find parent
        double spinWheel = 0;
        for (Individual individual : newIndividuals) {
            spinWheel += individual.getFitness();
            if (spinWheel >= rouletteWheelPosition) {
                return individual;
            }
        }

        return individuals[population.size() - 1];
    }

    /**
     * Apply crossover to population
     *
     * Crossover, more colloquially considered "mating", takes the population
     * and blends individuals to create new offspring. It is hoped that when two
     * individuals crossover that their offspring will have the strongest
     * qualities of each of the parents. Of course, it's possible that an
     * offspring will end up with the weakest qualities of each parent.
     *
     * This method considers both the GeneticAlgorithm instance's crossoverRate
     * and the elitismCount.
     *
     * The type of crossover we perform depends on the problem domain. We don't
     * want to create invalid solutions with crossover, so this method will need
     * to be changed for different types of problems.
     *
     * This particular crossover method selects random genes from each parent.
     *
     * @param population The population to apply crossover to
     * @return The new population
     */
    public Population crossoverPopulation(Population population) {
        // Create new population
        Population newPopulation = new Population(population.size());

        // Loop over current population by fitness
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual parent1 = population.getFittest(populationIndex);

            // Apply crossover to this individual?
            if (this.crossoverRate > Math.random() && populationIndex >= this.elitismCount) {
                // Initialize offspring
                Individual offspring = new Individual(parent1.getChromosomeLength());

                // Find second parent
                Individual parent2 = selectParentRank(population);

                // Loop over genome
                for (int geneIndex = 0; geneIndex < parent1.getChromosomeLength(); geneIndex++) {
                    // Use half of parent1's genes and half of parent2's genes
                    if (0.5 > Math.random()) {
                        offspring.setGene(geneIndex, parent1.getGene(geneIndex));
                    } else {
                        offspring.setGene(geneIndex, parent2.getGene(geneIndex));
                    }
                }

                // Add offspring to new population
                newPopulation.setIndividual(populationIndex, offspring);
            } else {
                // Add individual to new population without applying crossover
                newPopulation.setIndividual(populationIndex, parent1);
            }
        }

        return newPopulation;
    }

    public Population crossoverPopulationOnePoint(Population population) {
        // Create new population
        Population newPopulation = new Population(population.size());

        // Loop over current population by fitness
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual parent1 = population.getFittest(populationIndex);

            // Apply crossover to this individual?
            if (this.crossoverRate > Math.random() && populationIndex >= this.elitismCount) {
                // Initialize offspring
                Individual offspring = new Individual(parent1.getChromosomeLength());

                // Find second parent
                Individual parent2 = selectParentTournament(population);

                // Loop over genome
                if (0.5 > Math.random()) {
                    for (int i = 0; i < parent1.getChromosomeLength() / 2; i++) {
                        offspring.setGene(i, parent1.getGene(i));
                    }
                } else {
                    for (int i = 0; i < parent1.getChromosomeLength() / 2; i++) {
                        offspring.setGene(i, parent2.getGene(i));
                    }
                }

                if (0.5 > Math.random()) {
                    for (int i = parent1.getChromosomeLength() / 2; i < parent1.getChromosomeLength(); i++) {
                        offspring.setGene(i, parent1.getGene(i));
                    }
                } else {
                    for (int i = parent1.getChromosomeLength() / 2; i < parent1.getChromosomeLength(); i++) {
                        offspring.setGene(i, parent2.getGene(i));
                    }
                }

                // Add offspring to new population
                newPopulation.setIndividual(populationIndex, offspring);
            } else {
                // Add individual to new population without applying crossover
                newPopulation.setIndividual(populationIndex, parent1);
            }
        }

        return newPopulation;
    }

    public Population crossoverPopulationTwoPoint(Population population) {
        // Create new population
        Population newPopulation = new Population(population.size());

        // Loop over current population by fitness
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual parent1 = population.getFittest(populationIndex);

            // Apply crossover to this individual?
            if (this.crossoverRate > Math.random() && populationIndex >= this.elitismCount) {
                // Initialize offspring
                Individual offspring = new Individual(parent1.getChromosomeLength());

                // Find second parent
                Individual parent2 = selectParentTournament(population);

                // Loop over genome
                if (0.5 > Math.random()) {
                    for (int i = 0; i < parent1.getChromosomeLength() / 3; i++) {
                        offspring.setGene(i, parent1.getGene(i));
                    }
                } else {
                    for (int i = 0; i < parent1.getChromosomeLength() / 3; i++) {
                        offspring.setGene(i, parent2.getGene(i));
                    }
                }

                if (0.5 > Math.random()) {
                    for (int i = parent1.getChromosomeLength() / 3; i < 2 * parent1.getChromosomeLength() / 3; i++) {
                        offspring.setGene(i, parent1.getGene(i));
                    }
                } else {
                    for (int i = parent1.getChromosomeLength() / 3; i < 2 * parent1.getChromosomeLength() / 3; i++) {
                        offspring.setGene(i, parent2.getGene(i));
                    }
                }

                if (0.5 > Math.random()) {
                    for (int i = 2 * parent1.getChromosomeLength() / 3; i < parent1.getChromosomeLength(); i++) {
                        offspring.setGene(i, parent1.getGene(i));
                    }
                } else {
                    for (int i = 2 * parent1.getChromosomeLength() / 3; i < parent1.getChromosomeLength(); i++) {
                        offspring.setGene(i, parent2.getGene(i));
                    }
                }

                // Add offspring to new population
                newPopulation.setIndividual(populationIndex, offspring);
            } else {
                // Add individual to new population without applying crossover
                newPopulation.setIndividual(populationIndex, parent1);
            }
        }

        return newPopulation;
    }

    /**
     * Apply mutation to population
     *
     * Mutation affects individuals rather than the population. We look at each
     * individual in the population, and if they're lucky enough (or unlucky, as
     * it were), apply some randomness to their chromosome. Like crossover, the
     * type of mutation applied depends on the specific problem we're solving.
     * In this case, we simply randomly flip 0s to 1s and vice versa.
     *
     * This method will consider the GeneticAlgorithm instance's mutationRate
     * and elitismCount
     *
     * @param population The population to apply mutation to
     * @return The mutated population
     */
    public Population mutatePopulation(Population population) {
        // Initialize new population
        Population newPopulation = new Population(this.populationSize);

        // Loop over current population by fitness
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual individual = population.getFittest(populationIndex);

            // Loop over individual's genes
            for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
                // Skip mutation if this is an elite individual
                if (populationIndex > this.elitismCount) {
                    // Does this gene need mutation?
                    if (this.mutationRate > Math.random()) {
                        // Get new gene
                        int newGene = 1;
                        if (individual.getGene(geneIndex) == 1) {
                            newGene = 0;
                        }
                        // Mutate gene
                        individual.setGene(geneIndex, newGene);
                    }
                }
            }

            // Add individual to population
            newPopulation.setIndividual(populationIndex, individual);
        }

        // Return mutated population
        return newPopulation;
    }

}
