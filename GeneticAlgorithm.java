
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
/*
    Referensi algoritma genetik:
    https://github.com/Apress/genetic-algorithms-in-java-basics/tree/master/GA%20in%20Java/src/chapter2

    Referensi untuk web scraping test case:
    https://gist.github.com/korakot/5c8e21a5af63966d80a676af0ce15067
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
        int fitness = 0; //fitness awal-awal diinisialisasi 0

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
        int[][] move = { //arah pengecekan kotak
            {-1, -1},   //kiri atas
            {-1, 0},    //atas
            {-1, 1},    //kanan atas
            {0, -1},    //kiri
            {0, 0},     //kotaknya sendiri
            {0, 1},     //kanan
            {1, -1},    //kiri bawah
            {1, 0},     //bawah
            {1, 1}      //kanan bawah
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
    private boolean validate(int i, int j) {
        return i < n && i >= 0 && j < n && j >= 0;
    }

    /*  method untuk menghitung fitness dari populasi dan fitness setiap individual
        fitness populasi dihitung dengan cara menjumlahkan fitness dari semua individualnya  */
    public void evalPopulation(Population population) {
        int populationFitness = 0;   //fitness population diinisialisai dengan 0

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

    /*  method crossover random 
        untuk setiap gene dari individu baru dilakukan random untuk menentukan untuk
        gene tersebut, akan mengambil milik parent 1 atau parent 2. */
    public Population crossoverPopulationRandom(Population population) {
        Population newPopulation = new Population(population.size());   //inisialisasi populasi baru

        //loop untuk mengisi populasi dengan individual baru
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            //pilih parent pertama
            Individual parent1 = population.getFittest(populationIndex);    //individual terbaik ke-sekian
//            Individual parent1 = selectParentRoulette(population);          //menggunakan roulette wheel selection
//            Individual parent1 = selectParentRank(population);              //menggunakan rank selection
//            Individual parent1 = selectParentTournament(population);        //menggunakan tournament


            /*  lakukan crossover dan pemilihan secara elitism.
                pertama-tama lakukan pemilihan secara elitism terlebih dahulu sesuai jumlah yang ditentukan
                setelah itu lakukan crossover jika crossover ratenya berhasil. */
            if (this.crossoverRate > Math.random() && populationIndex >= this.elitismCount) {   //jika crossover rate berhasil
                Individual offspring = new Individual(parent1.getChromosomeLength());   //inisialisasi individu baru

                //pilih parent kedua
//                Individual parent2 = selectParentRoulette(population);        //menggunakan roulette wheel selection
                Individual parent2 = selectParentRank(population);            //menggunakan rank selection
//                Individual parent2 = selectParentTournament(population);      //menggunakan tournament

                //loop setiap gene dalam chromosome
                for (int geneIndex = 0; geneIndex < parent1.getChromosomeLength(); geneIndex++) {
                    //gunakan random untuk menentukan untuk gene diambil dari parent 1 atau parent 2
                    if (0.5 > Math.random()) {
                        offspring.setGene(geneIndex, parent1.getGene(geneIndex));   //ambil gene parent 1
                    } else {
                        offspring.setGene(geneIndex, parent2.getGene(geneIndex));   //ambil gene parent 2
                    }
                }

                //tambahkan individu baru ke populasi baru
                newPopulation.setIndividual(populationIndex, offspring);
            } else {
                //elitism, langsung tambahkan individu ke populasi baru
                newPopulation.setIndividual(populationIndex, parent1);
            }
        }

        return newPopulation;   //return populasi baru
    }

    //method crossover one-point
    public Population crossoverPopulationOnePoint(Population population) {
        Population newPopulation = new Population(population.size());   //inisialisasi populasi baru

        //loop untuk mengisi populasi dengan individual baru
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            //pilih parent pertama
            Individual parent1 = population.getFittest(populationIndex);    //individual terbaik ke-sekian
//            Individual parent1 = selectParentRoulette(population);          //menggunakan roulette wheel selection
//            Individual parent1 = selectParentRank(population);              //menggunakan rank selection
//            Individual parent1 = selectParentTournament(population);        //menggunakan tournament

            /* lakukan crossover dan pemilihan secara elitism.
                pertama-tama lakukan pemilihan secara elitism terlebih dahulu sesuai jumlah yang ditentukan
                setelah itu lakukan crossover jika crossover ratenya berhasil. */
            if (this.crossoverRate > Math.random() && populationIndex >= this.elitismCount) {   //jika crossover rate berhasil
                Individual offspring = new Individual(parent1.getChromosomeLength());   //inisialisasi individu baru

                //pilih parent kedua
//                Individual parent2 = selectParentRoulette(population);        //menggunakan roulette wheel selection
                Individual parent2 = selectParentRank(population);            //menggunakan rank selection
//                Individual parent2 = selectParentTournament(population);      //menggunakan tournament

                //random untuk gene pertama sampai gene ke-((panjang chromosome/2)-1)
                if (0.5 > Math.random()) {  //random unutuk memilih parent 1 atau parent 2
                    for (int i = 0; i < parent1.getChromosomeLength() / 2; i++) {  //pilih parent 1
                        offspring.setGene(i, parent1.getGene(i));   //copy gene parent 1 ke individu baru
                    }
                } else {
                    for (int i = 0; i < parent1.getChromosomeLength() / 2; i++) {   //pilih parent 2
                        offspring.setGene(i, parent2.getGene(i));   //copy gene parent 2 ke individu baru
                    }
                }

                //random untuk gene ke-(panjang chromosome/2) sampai gene terakhir
                if (0.5 > Math.random()) {  //random unutuk memilih parent 1 atau parent 2
                    for (int i = parent1.getChromosomeLength() / 2; i < parent1.getChromosomeLength(); i++) {   //pilih parent 1
                        offspring.setGene(i, parent1.getGene(i));   //copy gene parent 1 ke individu baru
                    }
                } else {
                    for (int i = parent1.getChromosomeLength() / 2; i < parent1.getChromosomeLength(); i++) {   //pilih parent 2
                        offspring.setGene(i, parent2.getGene(i));   //copy gene parent 2 ke individu baru
                    }
                }

                //tambahkan individu baru ke populasi baru
                newPopulation.setIndividual(populationIndex, offspring);
            } else {
                //elitism, langsung tambahkan individu ke populasi baru
                newPopulation.setIndividual(populationIndex, parent1);
            }
        }

        return newPopulation;   //return populasi baru
    }

    //method crossover two-point
    public Population crossoverPopulationTwoPoint(Population population) {
        Population newPopulation = new Population(population.size());    //inisialisasi populasi baru

        //loop untuk mengisi populasi dengan individual baru
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            //pilih parent pertama
            Individual parent1 = population.getFittest(populationIndex);    //individual terbaik ke-sekian
//              Individual parent1 = selectParentRoulette(population);          //menggunakan roulette wheel selection
//              Individual parent1 = selectParentRank(population);              //menggunakan rank selection
//              Individual parent1 = selectParentTournament(population);        //menggunakan tournament

            /* lakukan crossover dan pemilihan secara elitism.
                pertama-tama lakukan pemilihan secara elitism terlebih dahulu sesuai jumlah yang ditentukan
                setelah itu lakukan crossover jika crossover ratenya berhasil. */
            if (this.crossoverRate > Math.random() && populationIndex >= this.elitismCount) {   //jika crossover rate berhasil
                Individual offspring = new Individual(parent1.getChromosomeLength());   //inisialisasi individu baru

                //pilih parent kedua
//                Individual parent2 = selectParentRoulette(population);        //menggunakan roulette wheel selection
                Individual parent2 = selectParentRank(population);            //menggunakan rank selection
//                Individual parent2 = selectParentTournament(population);      //menggunakan tournament

                //random untuk gene pertama sampai gene ke-((panjang chromosome/3)-1)
                if (0.5 > Math.random()) {  //random unutuk memilih parent 1 atau parent 2
                    for (int i = 0; i < parent1.getChromosomeLength() / 3; i++) {   //pilih parent 1
                        offspring.setGene(i, parent1.getGene(i));   //copy gene parent 1 ke individu baru
                    }
                } else {
                    for (int i = 0; i < parent1.getChromosomeLength() / 3; i++) {   //pilih parent 2
                        offspring.setGene(i, parent2.getGene(i));   //copy gene parent 2 ke individu baru
                    }
                }

                //random untuk gene ke-(panjang chromosome/3) sampai gene ke-((2*panjang chromosome/3)-1)
                if (0.5 > Math.random()) {  //random unutuk memilih parent 1 atau parent 2
                    for (int i = parent1.getChromosomeLength() / 3; i < 2 * parent1.getChromosomeLength() / 3; i++) {   //pilih parent 1
                        offspring.setGene(i, parent1.getGene(i));   //copy gene parent 1 ke individu baru
                    }
                } else {
                    for (int i = parent1.getChromosomeLength() / 3; i < 2 * parent1.getChromosomeLength() / 3; i++) {   //pilih parent 2
                        offspring.setGene(i, parent2.getGene(i));   //copy gene parent 2 ke individu baru
                    }
                }

                //random untuk gene ke-(2*panjang chromosome/3) sampai gene terakhir
                if (0.5 > Math.random()) {  //random unutuk memilih parent 1 atau parent 2
                    for (int i = 2 * parent1.getChromosomeLength() / 3; i < parent1.getChromosomeLength(); i++) {   //pilih parent 1
                        offspring.setGene(i, parent1.getGene(i));   //copy gene parent 1 ke individu baru   
                    }
                } else {
                    for (int i = 2 * parent1.getChromosomeLength() / 3; i < parent1.getChromosomeLength(); i++) {   //pilih parent 2
                        offspring.setGene(i, parent2.getGene(i));   //copy gene parent 2 ke individu baru
                    }
                }

                //tambahkan individu baru ke populasi baru
                newPopulation.setIndividual(populationIndex, offspring);
            } else {
                //elitism, langsung tambahkan individu ke populasi baru
                newPopulation.setIndividual(populationIndex, parent1);
            }
        }

        return newPopulation;    //return populasi baru
    }

    //method untuk melakukan mutasi pada populasi
    public Population mutatePopulation(Population population) {
        //jika mutation rate berhasil lakukan mutasi
        if (this.mutationRate > Math.random()) {
            int index = (int) (Math.random() * populationSize);         //random individu yang akan di mutasi
            Individual individual = population.getIndividual(index);    //ambil individual

            int offset = (int) (Math.random() * (n * n)); //random gene yang dimutasi
            if (individual.getGene(offset) == 1) {
                individual.setGene(offset, 0);  //jika gene 1 ubah menjadi 0 
            } else {
                individual.setGene(offset, 1);  //jika gene 0 ubah menjadi 1
            }
        }

        return population;  //return populasi
    }
}
