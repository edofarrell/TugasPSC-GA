/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

/**
 *
 * @author
 */
public class Individual {

    private int[] chromosome;
    private double fitness = -1;

    /**
     * Initializes individual with specific chromosome
     *
     * @param chromosome The chromosome to give individual
     */
    public Individual(int[] chromosome) {
        // Create individual chromosome
        this.chromosome = chromosome;
    }

    /**
     * Initializes random individual.
     *
     * This constructor assumes that the chromosome is made entirely of 0s and
     * 1s, which may not always be the case, so make sure to modify as
     * necessary. This constructor also assumes that a "random" chromosome means
     * simply picking random zeroes and ones, which also may not be the case
     * (for instance, in a traveling salesman problem, this would be an invalid
     * solution).
     *
     * @param chromosomeLength The length of the individuals chromosome
     */
    public Individual(int chromosomeLength) {

        this.chromosome = new int[chromosomeLength];
        for (int gene = 0; gene < chromosomeLength; gene++) {
            if (0.5 < Math.random()) {
                this.setGene(gene, 1);
            } else {
                this.setGene(gene, 0);
            }

        }

    }

    /**
     * Gets individual's chromosome
     *
     * @return The individual's chromosome
     */
    public int[] getChromosome() {
        return this.chromosome;
    }

    /**
     * Gets individual's chromosome length
     *
     * @return The individual's chromosome length
     */
    public int getChromosomeLength() {
        return this.chromosome.length;
    }

    /**
     * Set gene at offset
     *
     * @param gene
     * @param offset
     * @return gene
     */
    public void setGene(int offset, int gene) {
        this.chromosome[offset] = gene;
    }

    /**
     * Get gene at offset
     *
     * @param offset
     * @return gene
     */
    public int getGene(int offset) {
        return this.chromosome[offset];
    }

    /**
     * Store individual's fitness
     *
     * @param fitness The individuals fitness
     */
    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    /**
     * Gets individual's fitness
     *
     * @return The individual's fitness
     */
    public double getFitness() {
        return this.fitness;
    }

    /**
     * Display the chromosome as a string.
     *
     * @return string representation of the chromosome
     */
    public String toString() {
        String output = "";
        for (int gene = 0; gene < this.chromosome.length; gene++) {
            output += this.chromosome[gene];
        }
        return output;
    }
}
