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

    string runDir = "10m_batchrun";

    fprintf(stdout, "Task %d checking in!\n", irank);

    string fullRunPath = "$WORK/alice/"+runDir;
    struct stat dirChecker;
    if(stat(fullRunPath.c_str(), &dirChecker)!=0){
        fprintf(stdout, "Task %d fullRunPath: %s\n", irank, fullRunPath.c_str());
        mkdir(fullRunPath.c_str(), 0777);
    }else{
        fprintf(stdout, "stat returned 0 for Task %d\n", irank);
    }

    string outputDir = "batch_";
    ostringstream ss;
    ss <<  irank;
    outputDir += ss.str();

    string fullOutPath = fullRunPath + "/" + outputDir;
    if(stat(fullOutPath.c_str(), &dirChecker)!=0){
        mkdir(fullOutPath.c_str(), 0777);
    }


}
