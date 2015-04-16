#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>

using namespace std;

int main (int argc, char* argv[]){

    int ierr = 0, nrank, irank;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    string input_dir_root, input_dir, input_file, output_dir, output_file;
    string run_file = "histo_gen.C";
    string aliroot_cmd, aliroot_source="ali";
    int num_processors, first_dir, last_dir, range, per_processor;

    if(argc != 8){
        printf("Incorrect number of arguments.  Correct usage is: \n");
        printf("gen_histos(inputDirRoot, inputFile, outputDir, outputFile, NumProcessors, firstDir, lastDir)\n");
        return 0;
    }else{

        /***************************************
         **** STUFF TO CHANGE BETWEEN RUNS *****
         ***************************************/
        input_dir_root = argv[1];

        input_file = argv[2];
        output_dir = argv[3];
        output_file = argv[4];

        num_processors = atoi(argv[5]);
        first_dir = atoi(argv[6]);
        last_dir = atoi(argv[7]);

        range = last_dir - first_dir + 1;
        per_processor = ceil((double)(range)/(double)(num_processors));
 
        /***************************************/

        //fprintf(stdout, "Task %d checking in!\n", irank);

        struct stat dirChecker;
        if(stat(output_dir.c_str(), &dirChecker)!=0){
            mkdir(output_dir.c_str(), 0777);
        }
        for(int i = first_dir+(irank*per_processor); i < first_dir+((irank+1)*per_processor); i++){
            if(i > last_dir)break;
            ostringstream ss;
            ss << i;
            input_dir = input_dir_root + ss.str();

            string aliroot_cmd = "/work/03093/deepat/alice/aliroot/vAN-20140818/build/bin/tgt_linuxx8664gcc/aliroot -b -l -q \'/work/03094/jtblair/utaustin/aliroot/"+run_file+"(\""+input_dir+"\", \""+input_file+"\",\""+output_dir+"\", \""+output_file+"\")\'";
            system(aliroot_source.c_str());
            system(aliroot_cmd.c_str());
        }
        ierr = MPI_Finalize();
        return ierr;
    }
