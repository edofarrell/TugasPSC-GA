
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

    private int n; //ukuran papan (nxn)
    /*
        array 2d untuk menyimpan papan permainan beserta angka pada setiap kotak.
        Jika kotak kosong maka akan diisi dengan -1.
    */
    private int[][] board;          
    
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

    /*
        Fitness sebuah individual dihitung dengan cara menjumlahkan setiap nilai dari:
            Untuk setiap kotak pada papan permainan lihat semua kotak di sekelilingnya 
            termasuk diagonal dan dirinya sendiri, kemudian hitung ada berapa kotak yang diisi. 
    
            Jika jumlahnya melebihi nilai yang berada pada kotak tersebut maka individual 
            tersebut tidak valid. Jadi fitnessnya dibuat 0.
    */
    //method menghitung fitness suatu individual
    public double calcFitness(Individual individual) {
        double fitness = 0; //fitness awal-awal diinisialisasi 0

        for (int i = 0; i < n; i++) {       //loop untuk setiap bari pada papan
            for (int j = 0; j < n; j++) {   //loop untuk setiap kolom pada papan
                if (board[i][j] >= 0) {     //jika kotak memiliki nilai 
                    int correctGenes = countNeighbour(i, j, individual);   //hitung kotak yang diisi di sekelilingnya
                    if (correctGenes <= board[i][j]) {  //jika jumlah masih valid
                        fitness += correctGenes;        //tambahkan ke fitness
                    } else {                            //jika jumlah melebihi nilai kotak yang bersangkutan
                        fitness = 0;                    //set fitness menjadi 0
                        individual.setFitness(fitness); //set fitness individual yang bersangkutan
                        return fitness;                 //return fitness
                    }
                }
            }
        }

        individual.setFitness(fitness); //set fitness individual yang bersangkutan
        return fitness;                 //return fitness
    }

    //method untuk menghitung banyak kotak yang diisi di sekeliling satu kotak
    private int countNeighbour(int i, int j, Individual individual) {
        int[][] move = {    //arah pengecekan kotak
            {-1, -1},       //kiri atas
            {-1, 0},        //atas
            {-1, 1},        //kanan atas
            {0, -1},        //kiri
            {0, 0},         //kotaknya sendiri
            {0, 1},         //kanan
            {1, -1},        //kiri bawah
            {1, 0},         //bawah
            {1, 1}          //kanan bawah
        };

        int count = 0;      //variabel untuk menghitung jumlah kotak yang diisi
        for (int k = 0; k < move.length; k++) { //loop untuk memeriksa semua kotak di sekeliling
            int newI = i + move[k][0];      //baris kotak yang inign diperiksa
            int newJ = j + move[k][1];      //kolom kotak yang ingin diperiksa
            if (validate(newI, newJ)) {     //cek apakah baris dan kolom masih berada dalam papan
                int offset = n * newI + newJ;   //index pada chromosome
                if (individual.getGene(offset) == 1) {  //jika kotak diisi
                    count++;                            //maka tambahkan count
                }
            }
        }

        return count;   //return hasil
    }

    //method untuk memeriksa apakah kotak pada baris ke-i, kolom ke-j masih valid (berada dalam papan)
    private boolean validate(int i, int j) { return i < n && i >= 0 && j < n && j >= 0; }

    /*  method untuk menghitung fitness dari populasi dan fitness setiap individual
        fitness populasi dihitung dengan cara menjumlahkan fitness dari semua individualnya  */
    public void evalPopulation(Population population) {
        double populationFitness = 0;   //fitness population diinisialisai dengan 0

        //loop untuk menghitung fitness dari setaip individual pada populasi
        for (Individual individual : population.getIndividuals()) {
            populationFitness += calcFitness(individual);   //hitung fitness dari individual dan tambahkan ke fitness populasi
        }

        population.setPopulationFitness(populationFitness); //set fitness population
    }

    /*  method untuk memeriksa apakah kondisi berhenti algoritma sudah terpenuhi atau belum.
        pada algoritma ini menggunakan banyak generasi sebagai kondisi berhenti.  */
    public boolean isTerminationConditionMet(int generation) {
        return generation >= this.numOfGeneration; //jika sudah mencapai target jumlah generasi return true
    }

    //method parent selection secara roulette wheel selection
    public Individual selectParentRoulette(Population population) {
        Individual individuals[] = population.getIndividuals(); //ambil semua individual dari populasi

        double populationFitness = population.getPopulationFitness();       //ambil fitness populasi
        double rouletteWheelPosition = Math.random() * populationFitness;   //putar roulette wheel

        //cari parent yang terpilih
        double currArea = 0;   //batas area individu saat ini
        for (Individual individual : individuals) {    //loop untuk setiap individual
            currArea += individual.getFitness();       //area individu di roulette wheel
            if (currArea >= rouletteWheelPosition) {   //jika hasil putaran berada dalam area individu
                return individual;                     //return individu yang bersangkutan
            }
        }
        return individuals[population.size() - 1];      //jika tidak ditemukan dalam loop, artinya yang terpilih adalah individu terakhir
    }

    //method parent selection secara rank selection
    public Individual selectParentRank(Population population) {
        Individual individuals[] = population.getIndividuals(); //ambil semua individual dari populasi
        
        Arrays.sort(individuals, new Comparator<Individual>() { //urutkan individual berdasarkan fitness dari besar ke kecil
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

        /* 
            hitung rank total dari semua individu. 
            menggunakan rumus baris dan deret: n/2 * (1+n)
            contoh: jika ada 4 individu maka hasilnya adalah 1+2+3+4.
        */
        double totalRank = (individuals.length / 2.0) * (1 + individuals.length);
        
        double area[] = new double[individuals.length];     //array untuk menyimpan area setiap individu di roulette wheel berdasarkan rank
        for (int i = individuals.length - 1; i >= 0; i--) { //loop untuk menghitung area setiap inidividu
            area[i] = 1.0 * (individuals.length - i) / totalRank;
        }

        double rouletteWheelPosition = Math.random();   //putar roulette wheel

        //cari parent yang terpilih
        double currArea = 0;    //batas area individu saat ini
        for (int i = 0; i < area.length; i++) {     //loop untuk setiap individual
            currArea += area[i];                    //area individu di roulette wheel
            if (currArea >= rouletteWheelPosition) {//jika hasil putaran berada dalam area individu
                return individuals[i];              //return individu yang bersangkutan
            }
        }

        return individuals[population.size() - 1];  //jika tidak ditemukan dalam loop, artinya yang terpilih adalah individu terakhir
    }

    //method parent selection secara tournament
    public Individual selectParentTournament(Population population) {
        Individual individuals[] = population.getIndividuals(); //ambil semua individual dari populasi
        
        Arrays.sort(individuals, new Comparator<Individual>() { //urutkan individual berdasarkan fitness dari besar ke kecil
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

        int memberPerGroup = 5;             //tentukan banyak member per group yang ingin ditandingkan
        double newPopulationFitness = 0;    //fitness populasi yang baru setelah dilakukan pertandingan

        //array untuk menyimpan individual yang memenangkan pertandingan
        Individual newIndividuals[] = new Individual[individuals.length / memberPerGroup]; 
        for (int i = 0; i < newIndividuals.length; i++) {   //loop sebanyak individual yang menang
            //ambil individual yang menang (dapat dilakukan seperti ini karena sudah diurutkan)
            newIndividuals[i] = individuals[i * memberPerGroup];    
            newPopulationFitness += newIndividuals[i].getFitness(); //tambahkan fitness inidividual yang menang ke fitness populasi
        }

        double rouletteWheelPosition = Math.random() * newPopulationFitness;    //putar roulette wheel

        //cari parent yang terpilih
        double currArea = 0;    //batas area individu saat ini
        for (Individual individual : newIndividuals) {  //loop untuk setiap individual
            currArea += individual.getFitness();        //area individu di roulette wheel
            if (currArea >= rouletteWheelPosition) {    //jika hasil putaran berada dalam area individu
                return individual;                      //return individu yang bersangkutan
            }
        }

        return individuals[population.size() - 1];  //jika tidak ditemukan dalam loop, artinya yang terpilih adalah individu terakhir
    }

    //method crossover random
    public Population crossoverPopulationRandom(Population population) {
        Population newPopulation = new Population(population.size());   //inisialisasi populasi baru

        //loop untuk 
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

    //method crossover one-point
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

    //method crossover two-point
    public Population crossoverPopulationTwoPoint(Population population) {
        // Create new population
        Population newPopulation = new Population(population.size());

        // Loop over current population by fitness
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual parent1 = population.getFittest(populationIndex);
//            Individual parent1 = selectParentRank(population);

            // Apply crossover to this individual?
            if (this.crossoverRate > Math.random() && populationIndex >= this.elitismCount) {
                // Initialize offspring
                Individual offspring = new Individual(parent1.getChromosomeLength());

                // Find second parent
                Individual parent2 = selectParentRank(population);

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
