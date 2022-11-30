
import java.util.Scanner;

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
        Scanner sc = new Scanner(System.in);
        int populationSize = 1000;
        double mutationRate = 0.001;
        double crossoverRate = 0.8;
        int elitismCount = 2;

        int maxFitness = 0;
//        int n = sc.nextInt();
        int n = 5;
        int[][] board = new int[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                board[i][j] = sc.nextInt();
                if (board[i][j] != -1) {
                    maxFitness += board[i][j];
                }
            }
        }
        int numOfGeneration = 1000;

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
            System.out.println("Best solution: " + population.getFittest(0).toString() + " " + population.getFittest(0).getFitness() + "/" + maxFitness);

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
