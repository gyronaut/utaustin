#include <string>
#include <iostream>
using namespace std;

void diagReader(string filename){

    FILE * file = fopen(filename.c_str(), "r");
    char line[256];
    char section[8];
    Float_t time=0, min=0, getevent_time=0, io_time=0, analysis_time=0, write_time=0; 
      

    while(fgets(&line, 256, file)){
//        fprintf(stdout, "%s\n", line);
        if(!sscanf(&line[0], "%s %*s %*s 0:%f:%f %*s", &section, &min, &time))continue;
        if(section[0]=='G'){
//            fprintf(stdout, "GETEVENT\n");
            getevent_time += time;
        }else if(section[0]=='A'){
//            fprintf(stdout, "ANALYSIS\n");
            analysis_time += time;
        }else if(section[0]=='I'){
//            fprintf(stdout, "IO\n");
            io_time += time;
        }else if(section[0]=='W'){
            fprintf(stdout, "%s WRITE time: %f\n", line, time);
            write_time += time + 60*min;
        }
    }

    fclose(file);
    
    fprintf(stdout, "Time for Reading: %f\n", io_time);
    fprintf(stdout, "Time for getting Event: %f\n", getevent_time);
    fprintf(stdout, "Time for analysis: %f\n", analysis_time);
    fprintf(stdout, "Time for writing: %f\n", write_time);
    fprintf(stdout, "Total Time for 100k Events: %f\n", (io_time + getevent_time + analysis_time + write_time));
    
}
