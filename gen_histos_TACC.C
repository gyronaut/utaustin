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
    string inputDir = "/work/03094/jtblair/alice/25082014_10m_batchrun/batch_";
    
    ostringstream ss;
    ss <<  irank;
    inputDir += ss.str();

    string inputFile = "25k_minbias.root";

    string runFile = "histo_gen.C";
    
    string outputDir = "/work/03094/jtblair/alice/25082014_10m_batchrun/histograms";
    string outputFile = "histo_"+ss.str()+".root";

    /***************************************/

    //fprintf(stdout, "Task %d checking in!\n", irank);

    struct stat dirChecker;
    if(stat(outputDir.c_str(), &dirChecker)!=0){
        mkdir(outputDir.c_str(), 0777);
    }
    string alirootCmd = "/work/03093/deepat/alice/aliroot/vAN-20140818/build/bin/tgt_linuxx8664gcc/aliroot -b -q \'./aliroot/"+runFile+"(\""+inputDir+"\", \""+inputFile+"\",\""+outputDir+"\", \""+outputFile+"\")\'";
    system(alirootCmd.c_str());

    ierr = MPI_Finalize();
    return ierr;
}
