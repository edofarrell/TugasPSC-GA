
/**
 *
 * @author
 *  Kelompok 15:
 *      Keannen Renaldo Halim - 6182001007
 *      Neil Christopher - 6182001010
 *      Edo Farrell Haryanto - 6182001025
 */
public class Individual {

    /*
        Papan berukuran nxn dimodelkan kedalam chromosome dengan bentuk array 1d dengan panjang nxn.
        Di mana cell di papan pada baris i dan kolom j ditempatkan pada array index ke-(i*n+j).
        
        Setiap gene menandakan apakah kotak yang bersangkutan di papan diisi atau tidak.
    
        Alel terdiri dari 0/1
        0 menandakan kotak yang bersangkutan tidak diisi.
        1 menandakan kotak yang bersangkutan diisi.
    */
    
    private int[] chromosome;       //chromosome dalam bentuk array 1d
    private double fitness = -1;    //fitness individual

    //constructor dengan input chromosome yang sudah jadi
    public Individual(int[] chromosome) {
        this.chromosome = chromosome;   //copy chromosome dari input
    }

    //constructor yang mengenerate chromosome baru secara random dengan panjang tertentu
    public Individual(int chromosomeLength) {
        this.chromosome = new int[chromosomeLength];            //inisialisasi chromosome (array)
        for (int gene = 0; gene < chromosomeLength; gene++) {   //loop untuk mengisi setiap gene dalam chromosome
            if (0.5 < Math.random()) {  //random isi gene 0/1
                this.setGene(gene, 1);  //isi gene 1
            } else {
                this.setGene(gene, 0);  //isi gene 0
            }
        }
    }

    //method getter chromosome
    public int[] getChromosome() { return this.chromosome; }

    //method getter panjang chromosome
    public int getChromosomeLength() { return this.chromosome.length; }

    //method setter gene chromosome pada index tertentu
    public void setGene(int offset, int gene) { this.chromosome[offset] = gene; }

    //method getter gene chromose pada index tertentu
    public int getGene(int offset) { return this.chromosome[offset]; }

    //method setter fitness individual
    public void setFitness(double fitness) { this.fitness = fitness; }

    //method getter fitness individual
    public double getFitness() { return this.fitness; }

    //method untuk mengubah indivual menjadi bentuk string
    @Override
    public String toString() {
        String output = "";
        for (int gene = 0; gene < this.chromosome.length; gene++) {
            output += this.chromosome[gene];
        }
        return output;
    }
}
