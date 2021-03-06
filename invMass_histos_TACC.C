#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

int main (int argc, char* argv[]){

    int ierr = 0, nrank, irank;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    /***************************************
     **** STUFF TO CHANGE BETWEEN RUNS *****
     ***************************************/
    string inputDir = "/scratch/03094/jtblair/alice/20140912_100m_batchrun_1/batch_";
    
    ostringstream ss;
    ss <<  irank;
    inputDir += ss.str();

    string inputFile = "25k_minbias.root";

    string runFile = "k_invariant_mass.C";
    
    string outputDir = "/scratch/03094/jtblair/alice/20140912_100m_batchrun_1/invMass_histograms";
    string outputFile = "invMass_histo_"+ss.str()+".root";

    /***************************************/

    //fprintf(stdout, "Task %d checking in!\n", irank);

    struct stat dirChecker;
    if(stat(outputDir.c_str(), &dirChecker)!=0){
        mkdir(outputDir.c_str(), 0777);
    }
    string alirootCmd = "/work/03093/deepat/alice/aliroot/vAN-20140818/build/bin/tgt_linuxx8664gcc/aliroot -b -l -q \'/work/03094/jtblair/utaustin/aliroot/"+runFile+"(\""+inputDir+"\", \""+inputFile+"\",\""+outputDir+"\", \""+outputFile+"\")\'";
    system(alirootCmd.c_str());

    ierr = MPI_Finalize();
    return ierr;
}
