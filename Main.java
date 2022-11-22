/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

/**
 *
 * @author
 */
public class Main {

    public static void main(String[] args) {
        int populationSize = 100;
        double mutationRate = 0.001;
        double crossoverRate = 0.95;
        int elitismCount = 2;

        int n = 5;
        int[][] board = {
            {1, -1, 0, -1, -1},
            {-1, -1, -1, 3, -1},
            {-1, 5, -1, -1, 5},
            {-1, -1, -1, -1, -1},
            {-1, 4, 5, -1, 2}
        };
        int numOfGeneration = 100;

        int chromosomeLength = n * n;

        // Create GA object
        GeneticAlgorithm ga = new GeneticAlgorithm(board, numOfGeneration, populationSize, mutationRate, crossoverRate, elitismCount);

        // Initialize population
        Population population = ga.initPopulation(chromosomeLength);

        // Evaluate population
        ga.evalPopulation(population);

        // Keep track of current generation
        int generation = 0;

        /**
         * Start the evolution loop
         *
         * Every genetic algorithm problem has different criteria for finishing.
         * In this case, we know what a perfect solution looks like (we don't
         * always!), so our isTerminationConditionMet method is very
         * straightforward: if there's a member of the population whose
         * chromosome is all ones, we're done!
         */
//        while (ga.isTerminationConditionMet(population) == false) {
//            // Print fittest individual from population
//            System.out.println("Best solution: " + population.getFittest(0).toString());
//
//            // Apply crossover
//            population = ga.crossoverPopulation(population);
//
//            // Apply mutation
//            population = ga.mutatePopulation(population);
//
//            // Evaluate population
//            ga.evalPopulation(population);
//
//            // Increment the current generation
//            generation++;
//        }
        while (ga.isTerminationConditionMet(generation) == false) {
            // Print fittest individual from population
            System.out.println("Best solution: " + population.getFittest(0).toString() + " " + population.getFittest(0).getFitness());

//            System.out.println("Sorted\n");
//            for (int i = 0; i < populationSize; i++) {
//                System.out.println(population.getIndividual(i).toString() + " " + population.getIndividual(i).getFitness());
//            }
//            System.out.println("");
            // Apply crossover
            population = ga.crossoverPopulation(population);

            // Apply mutation
            population = ga.mutatePopulation(population);

            // Evaluate population
            ga.evalPopulation(population);

            // Increment the current generation
            generation++;
        }

        /**
         * We're out of the loop now, which means we have a perfect solution on
         * our hands. Let's print it out to confirm that it is actually all
         * ones, as promised.
         */
        System.out.println("Found solution in " + generation + " generations");
        System.out.println("Best solution: " + population.getFittest(0).toString());
    }
}
