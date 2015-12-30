#include <cstdio>
#include <cstdlib>
#include <cmath>

int main(int argc,char** argv){
    if(argc!=2){
        printf("Usage: %s (m)",argv[0]);
    }
    double output=0;
    double inputM=atof(argv[1]);
    for(double k=1;k<10000000;k++){
        output+=1/k/(k+inputM);
    }
    printf("%2.10lf\n",output);
    return 0;
}
