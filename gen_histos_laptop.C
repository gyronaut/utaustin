#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>

using namespace std;

int main (int argc, char* argv[]){

    string input_dir_root, input_dir, input_file, output_dir, output_file_root, output_file;
    string run_file = "/Users/jtblair/utaustin/aliroot/K0_histo_gen.C";
    string aliroot_cmd, aliroot_source="ali";
    int num_processors, first_dir, last_dir, range, per_processor;

    if(argc != 7){
        printf("Incorrect number of arguments.  Correct usage is: \n");
        printf("gen_histos(input_dir_root, input_file, output_dir, output_file_root, first_dir, last_dir)\n");
        return 1;
    }else{

        input_dir_root.assign(argv[1], find(argv[1], argv[1]+255, '\0'));

        input_file.assign(argv[2], find(argv[2], argv[2]+255, '\0'));
        output_dir.assign(argv[3], find(argv[3], argv[3]+255, '\0'));
        output_file_root.assign(argv[4], find(argv[4], argv[4]+255, '\0'));

        first_dir = atoi(argv[5]);
        last_dir = atoi(argv[6]);

        range = last_dir - first_dir + 1;
        if(range <= 0){
            fprintf(stderr, "Error: first_dir must be less than or equal to last_dir number");
        }

        per_processor = ceil((double)(range)/(double)(num_processors));

        struct stat dirChecker;
        if(stat(output_dir.c_str(), &dirChecker)!=0){
            mkdir(output_dir.c_str(), 0777);
        }         

//        system(aliroot_source.c_str());

        for(int i = first_dir; i < last_dir; i++){
            ostringstream ss;
            ss << i;
            input_dir = input_dir_root + ss.str();

            output_file = output_file_root + ss.str()+".root";

            string aliroot_cmd = "aliroot -b -l -q \'"+run_file+"(\""+input_dir+"\", \""+input_file+"\",\""+output_dir+"\", \""+output_file+"\")\'";
            system(aliroot_cmd.c_str());
        }
    }
    return 0;
}
