
import java.util.Scanner;

/**
 *
 * @author 
 *  Kelompok 15: 
 *      Keannen Renaldo Halim - 6182001007 
 *      Neil Christopher - 6182001010 
 *      Edo Farrell Haryanto - 6182001025
 *
 * Sumber:
 *      https://github.com/Apress/genetic-algorithms-in-java-basics/tree/master/GA%20in%20Java/src/chapter2
 */

public class Main {

    public static void main(String[] args) {
        int populationSize = 1000;      //besar populasi
        double mutationRate = 0.001;    //probabilitas terjadi mutasi
        double crossoverRate = 0.8;     //probabilitas crossover berhasil
        int elitismCount = 2;           //jumlah individu yang akan dipilih secara elitism
        int numOfGeneration = 2200;     //banyak generasi

        Scanner sc = new Scanner(System.in);    //objek scanner
//        int n = sc.nextInt();
        int n = 5;                              //besar papan (nxn)
        int chromosomeLength = n * n;           //panjang kromosom = besar papan

        int buang = sc.nextInt();               //varibel untuk membuang 1 angka pertma dari input karena input kelebihan 1 angka yang tidak digunakan

        int maxFitness = 0;                     //variabel untuk menghitung  nilai fitness maximum yang bisa didapatkan
        int[][] board = new int[n][n];          //array 2d untuk menyimpan papan permainan (setiap cell diisi dengan angka, jika cell tersebut kosong maka diisi dengan -1)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                board[i][j] = sc.nextInt();     //mengisi papan permainan
                if (board[i][j] != -1) {
                    maxFitness += board[i][j];  //menghitung nilai fitness maksimum
                }
            }
        }

        // Membuat object GeneticAlgorithm
        GeneticAlgorithm ga = new GeneticAlgorithm(board, numOfGeneration, populationSize, mutationRate, crossoverRate, elitismCount);
        Population population = ga.initPopulation(chromosomeLength);    //inisialisasi populasi
        ga.evalPopulation(population);  //evaluasi populasi
        int generation = 0;             //varibel untuk menyimpan nomor generasi saat ini

        //loop selama syarat terminasi belum terpenuhi (jumlah generasi)
        while (ga.isTerminationConditionMet(generation) == false) {
            // print individual terbaik dari populasi
            System.out.println("Best solution: " + population.getFittest(0).toString() + " " + population.getFittest(0).getFitness() + "/" + maxFitness);

            population = ga.crossoverPopulationTwoPoint(population);    //crossover populasi
            population = ga.mutatePopulation(population);   //mutasi populasi
            ga.evalPopulation(population);                  //evaluasi populasi
            generation++;   //lanjut ke generasi berikutnya
        }

        //print hasil
        System.out.println("Found solution in " + generation + " generations");
        System.out.println("Best solution: " + population.getFittest(0).toString());
    }
}
