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
    string runDir = "20150302_100m_batchrun_1";
    
    string outputDir = "batch_";
    ostringstream ss;
    ss <<  irank;
    outputDir += ss.str();

    int nEvts = 25000;
    ostringstream ss2;
    ss2 << nEvts;
    string numEvents = ss2.str();
    
    string runFile = "minbias_batch_run.C";

    string outFile = "25k_minbias.root";
    /***************************************/

    //fprintf(stdout, "Task %d checking in!\n", irank);

    string fullRunPath = "/scratch/03094/jtblair/alice/"+runDir;
    struct stat dirChecker;
    if(stat(fullRunPath.c_str(), &dirChecker)!=0){
//      fprintf(stdout, "Task %d fullRunPath: %s\n", irank, fullRunPath.c_str());
        mkdir(fullRunPath.c_str(), 0777);
    }else{
//      fprintf(stdout, "stat returned 0 for Task %d\n", irank);
    }

    string fullOutPath = fullRunPath + "/" + outputDir+"/";
    if(stat(fullOutPath.c_str(), &dirChecker)!=0){
        mkdir(fullOutPath.c_str(), 0777);
    }

    string alirootCmd = "/work/03093/deepat/alice/aliroot/vAN-20140818/build/bin/tgt_linuxx8664gcc/aliroot -b -q \'/work/03094/jtblair/utaustin/aliroot/"+runFile+"("+numEvents+", \""+fullOutPath+outFile+"\", \"FASTGEN\")\' > /dev/null";

//  fprintf(stdout, "%s", alirootCmd.c_str());
    system(alirootCmd.c_str());

    ierr = MPI_Finalize();
    return ierr;
}
