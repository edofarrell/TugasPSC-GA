
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;
import java.io.FileReader;
import com.opencsv.CSVWriter;
import java.io.FileNotFoundException;

/**
 *
 * @author 
 *  Kelompok 15: 
 *      Keannen Renaldo Halim   - 6182001007
 *      Neil Christopher        - 6182001010 
 *      Edo Farrell Haryanto    - 6182001025
 */
public class Main {

    public static void main(String[] args) throws FileNotFoundException {
        int numOfGeneration = 1000;     // banyak generasi
        int populationSize = 1000;      // besar populasi
        double mutationRate = 0.001;    // probabilitas terjadi mutasi
        double crossoverRate = 0.8;     // probabilitas crossover berhasil
        int elitismCount = 10;          // jumlah individu yang akan dipilih secara elitism

        File folder = new File("./selenium/"+"test case"+"/");  // ambil folder yang berisi file input
        File[] input = folder.listFiles();                      // ambil semua file input

        for (int in = 0; in < input.length; in++) {     // kerjakan setiap test case
            String tc = "";                      // untuk menampung test case dalam bentuk String
            Scanner sc = new Scanner(input[in]); // objek scanner

            // int n = sc.nextInt();
            int n = 5; // besar papan (nxn)
            int chromosomeLength = n * n; // panjang kromosom = besar papan

            // varibel untuk membuang 1 angka pertma dari input karena 
            // input kelebihan 1 angka(dari web scrapping) yang tidak digunakan
            int buang = sc.nextInt(); 

            int maxFitness = 0;             // variabel untuk menghitung nilai fitness maximum yang bisa didapatkan
            int[][] board = new int[n][n];  // array 2d untuk menyimpan papan permainan (setiap cell diisi dengan angka,
                                            // jika cell tersebut kosong maka diisi dengan -1)
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    board[i][j] = sc.nextInt();     // mengisi papan permainan
                    tc += board[i][j] + " ";
                    if (board[i][j] != -1) {
                        maxFitness += board[i][j];  // menghitung nilai fitness maksimum
                    }
                }
                tc += "\n";
            }

            // Membuat object GeneticAlgorithm
            GeneticAlgorithm ga = new GeneticAlgorithm(board, numOfGeneration, populationSize, mutationRate, crossoverRate, elitismCount);
            Population population = ga.initPopulation(chromosomeLength); // inisialisasi populasi
            ga.evalPopulation(population);  // evaluasi populasi
            int generation = 0;             // varibel untuk menyimpan nomor generasi saat ini

            // loop selama syarat terminasi belum terpenuhi (jumlah generasi)
            while (ga.isTerminationConditionMet(generation) == false) {
                // print individual terbaik dari populasi
                System.out.println("Best solution: " + population.getFittest(0).toString() + ", Fitness: "
                        + population.getFittest(0).getFitness() + "/" + maxFitness);

                population = ga.crossoverPopulationTwoPoint(population); // crossover populasi
                population = ga.mutatePopulation(population); // mutasi populasi
                ga.evalPopulation(population);  // evaluasi populasi
                generation++;                   // lanjut ke generasi berikutnya
            }

            // print hasil
            System.out.println("Found solution in " + generation + " generations");
            System.out.println("Best solution: " + population.getFittest(0).toString());

            //data untuk menyimpan ke file
            String filePath = "./result.csv";                   // path file output
            String PopulationSize = populationSize + "";        // besar populasi
            String MutationRate = String.format("%.4f", mutationRate);  // probabilitas terjadi mutasi
            String CrossoverRate = crossoverRate + "";          // probabilitas crossover berhasil
            String NumberOfGenerations = numOfGeneration + "";  // banyak generasi
            String NumberOfElitism = elitismCount + "";         // jumlah individu yang akan dipilih secara elitism
            String BestChromosome = population.getFittest(0).toString(); // individu terbaik
            String FitnessChromosome = population.getFittest(0).getFitness() + ""; // fitness individu terbaik
            String Maxfitness = maxFitness + "";    //nilai fitness maksimum untuk test case
            String TestCase = tc;                   //test case
            writeDataLineByLine(    //tulis ke file output
                    filePath,
                    PopulationSize,
                    MutationRate,
                    CrossoverRate,
                    NumberOfGenerations,
                    NumberOfElitism,
                    BestChromosome,
                    FitnessChromosome,
                    Maxfitness,
                    TestCase
            );
        }
    }

    //method untuk menulis hasil eksperimen pada file csv
    public static void writeDataLineByLine(
            String filePath,
            String PopulationSize,
            String MutationRate,
            String CrossoverRate,
            String NumberOfGenerations,
            String NumberOfElitism,
            String BestChromosome,
            String FitnessChromosome,
            String Maxfitness,
            String TestCase
    ) {
        // buat file untuk menyimpan hasil eksperimen
        File file = new File(filePath);
        try {

            FileWriter outputfile = new FileWriter(file, true); // object filewriter untuk menulis ke file output

            CSVWriter writer = new CSVWriter(outputfile);       // object csvwriter untuk menulis ke file output dalam format csv

            // penambahan header pada csv
            BufferedReader br = new BufferedReader(new FileReader(file));   //object bufferedreader
            if (br.readLine() == null) {        //jika file belum pernah diisi
                String[] header = {             //masukkan header
                    "Population Size",
                    "Mutation Rate",
                    "Crossover Rate",
                    "Number of Generations",
                    "Number of Elitism",
                    "Best Chromosome",
                    "Fitness Chromosome",
                    "Max fitness",
                    "Test Case"
                };
                writer.writeNext(header);   //masukkan header
            }

            //masukkan data eksperimen pada csv
            String[] data1 = {
                PopulationSize,
                MutationRate,
                CrossoverRate,
                NumberOfGenerations,
                NumberOfElitism,
                BestChromosome,
                FitnessChromosome,
                Maxfitness,
                TestCase
            };
            writer.writeNext(data1);    //tulis hasil eksperimen pada file csv

            writer.close(); //close csvwriter
        } catch (IOException e) {
        }
    }

}
