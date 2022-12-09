
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

public class Population {

    private Individual population[];        //array 1d untuk menyimpan semua individual di populasi
    private double populationFitness = -1;  //variabel untuk menyimpan fitness populasi

    //constructor yang menginisialisasi populasi
    public Population(int populationSize) {
        this.population = new Individual[populationSize];  //inisialisasi array populasi
    }

    //constrcutor yang menginisialisasi populasi sekaligus setiap individualnya
    public Population(int populationSize, int chromosomeLength) {
        this.population = new Individual[populationSize];   //inisialisasi array populasi

        //loop untuk inisialisasi semua individual dalam populasi
        for (int individualCount = 0; individualCount < populationSize; individualCount++) {
            Individual individual = new Individual(chromosomeLength);   //inisialisasi individual
            this.population[individualCount] = individual;  //masukkan individual ke dalam populasi
        }
    }

    //method getter populasi
    public Individual[] getIndividuals() { return this.population; }

    //method getter untuk mengembalikan individual terbaik di populasi
    public Individual getFittest(int offset) {
        //sort populasi berdasarkan fitness individual dari besar ke kecil
        Arrays.sort(this.population, new Comparator<Individual>() {     
            @Override
            public int compare(Individual o1, Individual o2) {  //method compare
                if (o1.getFitness() > o2.getFitness()) {        //dari besar ke kecil
                    return -1;
                } else if (o1.getFitness() < o2.getFitness()) {
                    return 1;
                }
                return 0;
            }
        });

        return this.population[offset]; //return the individual dengan fitness terbaik
    }

    //method setter fitness population
    public void setPopulationFitness(double fitness) { this.populationFitness = fitness; }

    //method getter fitness population
    public double getPopulationFitness() { return this.populationFitness; }

    //method getter ukuran population
    public int size() { return this.population.length; }

    //method setter individual pada index tertentu
    public void setIndividual(int offset, Individual individual) { population[offset] = individual; }

    //method getter inidividual pada index tertentu
    public Individual getIndividual(int offset) { return population[offset]; }
}
