/* Include directives */
#include "ffpopsim_highd.h"
#include "multi_population.h"
#include <cstring>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cmath>
#include <boost/foreach.hpp>
#include <iostream>
#include <fstream>
#include <boost/format.hpp>

//#include <boost/filesystem.hpp>
#define HIGHD_BADARG -1354341
#define NOTHING 1e-10

/* Be verbose? */
#define HIGHD_VERBOSE 0


//#define MUTATION_RATE 1e-5
//#define MIGRATION_RATE 0
#define NUMBER_OF_TAITS 2
#define _L 100
#define GENERATIONS_Eq  150000
#define SHIFT_BEGIN 20000
#define SHIFT_END 90000
//#define SELECTIVE 1
//#define SHIFT_RATE 5000
#define GENERATIONS  10000
#define Av_Num 0
//#define SAMPLE_SIZE 1500


int main(int argc, char  *argv[]){

    //
    // ############################################
    ////
    //// Read the input parameters of the simulation
    ////
    //// ############################################
    //int SELECTIVE = atof(argv[1]); //MIGRATION_RATE = atof(argv[1]);

    int SEQ_LEN = atoi(argv[1]);
    int POPULATION_SIZE = atoi(argv[2]);
    int SAMPLE_NUM = atoi(argv[3]);
    int SAMPLE_FREQ = atoi(argv[4]);
    int SAMPLE_VOL = atoi(argv[5]);
    double MU = atof(argv[6]);
    std::string base_name = argv[7];

    std::cout << "SEQ_LEN = "  <<  SEQ_LEN   << std::endl;
    std::cout << "POPULATION_SIZE = " << POPULATION_SIZE << std::endl;
    std::cout << "SAMPLE_NUM = " << SAMPLE_NUM << std::endl;
    std::cout << "SAMPLE_FREQ = " << SAMPLE_FREQ << std::endl;
    std::cout << "SAMPLE_VOL = " << SAMPLE_VOL << std::endl;
    std::cout << "MU = " << MU << std::endl;

    //const char * fileDir =  argv[3];
    //int locations = atoi(argv[4]);
    //double MUTATION_RATE = atof(argv[5]);
    //double MIGRATION_RATE = atof(argv[6]);
    //// int _L = atoi(argv[5]);

    //int L = SEQ_LEN;
    //int N = POPULATION_SIZE;
    //int Ns = SAMPLE_NUM;
    //int Ts = SAMPLE_FREQ;

    //boost::format fmter("FFpopSim_L%1%_N%2%_Ns%3%_Ts%4%_Nv%5%_Mu%6%");

    haploid_highd pop(SEQ_LEN);
    pop.set_mutation_rate(MU);

    // genealogy
    vector <int> gen_loci;
    gen_loci.push_back(100);
    pop.track_locus_genealogy(gen_loci);

    pop.set_wildtype(POPULATION_SIZE);		// start with a population of the right size

    std::cout << "EQUILIBRATING POPULATION..." << std::endl;
    pop.set_mutation_rate(1e-2);
    pop.evolve(5*POPULATION_SIZE);
    // every site mutates ~5 times -> equilibrium
    pop.set_mutation_rate(MU);
    pop.evolve(5*POPULATION_SIZE);

    std::cout << "GENERATING TIME-RESOLVED SAMPLING..." << std::endl;
    for (int n=0; n<SAMPLE_NUM; n++){
        std::cout << n << " out of " << SAMPLE_NUM << std::endl;
        pop.evolve(SAMPLE_FREQ);
        pop.genealogy.trees[0].check_tree_integrity();
        pop.tree_sample=SAMPLE_VOL;
        pop.evolve(1); // take sample
        pop.genealogy.trees[0].check_tree_integrity();
        pop.tree_sample=0;
    }

    auto ts = pop.genealogy.trees;
    auto tree = ts.at(0);
    rooted_tree new_tree;
    tree.construct_subtree(tree.sampled_leafs, new_tree);

    /// save results
    //fmter % SEQ_LEN % POPULATION_SIZE % SAMPLE_NUM % SAMPLE_FREQ % SAMPLE_VOL % MU;
    //auto base_name = fmter.str();

    std::ofstream ofs;
    ofs.open(base_name + ".nwk");
    ofs << new_tree.print_newick() << std::endl;
    ofs.close();

    ofs.open(base_name + ".bin.fasta");
    ofs << tree.print_sequences() << std::endl;
    ofs.close();

    return 0;

}
