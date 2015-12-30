#include "gf_complete.h"
#include "NativeStair.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int main(int argc,char** argv){
  int i;
  int* globalParityVector=(int*)calloc(sizeof(int),2);
  globalParityVector[0]=1;
  globalParityVector[1]=2;
  nativeInit(9,6,3,5,1,globalParityVector,2);
  bool* erasedLoc=(bool*)calloc(sizeof(bool),15);
  char** input=(char**)calloc(sizeof(char*),15);
  char** output=(char**)calloc(sizeof(char*),15);
  int* readLayout=(int*)calloc(sizeof(int),15);
  int inputFd=open("input",O_RDONLY);
  for(i=0;i<15;i++){
    input[i]=(char*)calloc(sizeof(char),1024);
    if(i<9){
      pread(inputFd,input[i],1024,i*1024);
    }
    //if(i==0||i==6||i==7||i==8||i==9) output[i]=(char*)calloc(sizeof(char),1024);
    output[i]=(char*)calloc(sizeof(char),1024);
    erasedLoc[i]=false;
    readLayout[i]=0;
  }
  close(inputFd);
  erasedLoc[6]=true;
  //erasedLoc[7]=true;
  char** encodedBufs=(char**)calloc(sizeof(char*),15);
  for(i=0;i<6;i++){
    encodedBufs[i]=output[i];
  }
  for(i=0;i<9;i++){
    encodedBufs[i+6]=input[i];
  }
  nativeLocationsToReadForDecode(erasedLoc,readLayout,15);
  nativeEncodeBulk(input,9,output,6,1024);
  for(i=0;i<15;i++){
    printf("readLayout[%d]=%d\n",i,readLayout[i]);
  }
  char** recInput=(char**)calloc(sizeof(char*),15);
  char** recoveredOutput=(char**)calloc(sizeof(char*),15);

  // TODO: make the following piece of code work!!
  for(i=0;i<15;i++){
    if(readLayout[i]==1){
      recInput[i]=encodedBufs[i];
    }else if(erasedLoc[i]){
      recInput[i]=(char*)calloc(sizeof(char),1024);
    }else{
      //recInput[i]=(char*)calloc(sizeof(char),1024);
      recInput[i]=NULL;
    }
  }
  for(i=0;i<1;i++){
    recoveredOutput[i]=(char*)calloc(sizeof(char),1024);
  }
  nativeDecodeBulk(recInput,4,recoveredOutput,1,erasedLoc,1,1024);
  return 0;
}

