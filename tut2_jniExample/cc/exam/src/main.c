#include<stdio.h>
#include<stdlib.h>

#define MAX_LEN 100
#define SKIP_STATE 0
#define SAVE_STATE 1

int capitalizeStr(char* str,int len){
  int i;
  for(i=0;i<len;i++){
    if(str[i]>=97&&str[i]<=122){
      str[i]-=32;
    }
  }
  return 0;
}

int main(int argc,char** argv){
  int i;
  char* input=(char*)calloc(sizeof(char),MAX_LEN+1);
  char* output=(char*)calloc(sizeof(char),MAX_LEN+1);
  int inputIndex=0;
  int outputIndex=0;
  int currentState=0;
  fgets(input,MAX_LEN,stdin);
  printf("[pre-processed]: %s",input);
  while(1){
    if(input[inputIndex]=='\n'){
      break;
    }
    if(currentState==SKIP_STATE){
      if(input[inputIndex]==' '){
        inputIndex++;
        continue;
      }else{
        output[outputIndex]=input[inputIndex];
        inputIndex++;
        outputIndex++;
        currentState=SAVE_STATE;
      }
    }else{
      // save
      if(input[inputIndex]==' '){
        output[outputIndex]=input[inputIndex];
        inputIndex++;
        outputIndex++;
        currentState=SKIP_STATE;
      }else{
        output[outputIndex]=input[inputIndex];
        inputIndex++;
        outputIndex++;
      }
    }
  }
  if(output[outputIndex-1]==' '){
    output[outputIndex-1]='\0';
  }
  printf("ouput= %s\n",output);
  int prevSpaceLoc=-1;
  int maxSubStrBegining=-1;
  int maxSubStrEnding=-1;
  int maxSubStrLen=-1;
  for(i=0;i<outputIndex;i++){
    if(output[i]==' '){
      if(i-prevSpaceLoc-1>maxSubStrLen){
        printf("i=%d prev=%d maxLen=%d\n",i,prevSpaceLoc,maxSubStrLen);
        maxSubStrBegining=prevSpaceLoc+1;
        maxSubStrEnding=i-1;
        maxSubStrLen=i-prevSpaceLoc-1;
        printf(" i=%d prev=%d\n",i,prevSpaceLoc);
      }
      prevSpaceLoc=i;
    }if(i==outputIndex-1){
      if(i-prevSpaceLoc>maxSubStrLen){
        printf("i=%d prev=%d maxLen=%d\n",i,prevSpaceLoc,maxSubStrLen);
        maxSubStrBegining=prevSpaceLoc+1;
        maxSubStrEnding=i;
        maxSubStrLen=i-prevSpaceLoc;
        printf(" i=%d prev=%d\n",i,prevSpaceLoc);
      }
    }
  }
  char* final=(char*)calloc(sizeof(char),outputIndex+3);
  printf("longest sub-string: begins %d end %d\n",maxSubStrBegining,maxSubStrEnding);
  strncpy(final,output,maxSubStrBegining);
  final[maxSubStrBegining]='"';
  printf("[post-processing]: %s\n",final);
  strncpy(final+maxSubStrBegining+1,output+maxSubStrBegining,maxSubStrLen);
  final[maxSubStrEnding+2]='"';
  printf("[post-processing]: %s\n",final);
  strncpy(final+maxSubStrEnding+3,output+maxSubStrEnding+1,outputIndex-1-maxSubStrEnding);
  capitalizeStr(final+maxSubStrBegining+1,maxSubStrLen);
  printf("[post-processing]: %s\n",final);
  return 0;
}


