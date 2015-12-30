#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LEN 100

int main(){
  char current[MAX_LEN];
  int currentIdx=0;
  int inStringFlag=0;
  int count=0;
  int inputIdx=0;
  char input[MAX_LEN];
  int valid;
  int i;
  for(i=0;i<MAX_LEN;i++){
    input[i]=0;
  }
  fgets(input,MAX_LEN,stdin);
  for(i=0;i<MAX_LEN;i++){
    if(input[i]=='\n'){
      input[i]=0;
    }
  }
  //printf("%s\n",input);
  while(inputIdx<MAX_LEN) {
    if((input[inputIdx]==' ')||(input[inputIdx]==0)){
      current[currentIdx]=0;
      //printf("%s\n",current);
      if(inStringFlag==0){
        // do nothing!
        ;
      }else{
        current[currentIdx]=0;
        // judge whether current is a identifier
        if((current[0]<='9')&&(current[0]>='0')){
          ;
        }else if((currentIdx==2)&&
            (
             (strncmp(current,"if",2)==0)||
             (strncmp(current,"do",2)==0)
             )
            ){
          // not valid, skip
          ;
        }else if((currentIdx==3)&&
            (
             (strncmp(current,"int",3)==0)||
             (strncmp(current,"for",3)==0)
             )
            ){
          // not valid, skip
          ;
        }else if((currentIdx==4)&&
            (
             (strncmp(current,"auto",4)==0)||
             (strncmp(current,"char",4)==0)||
             (strncmp(current,"case",4)==0)||
             (strncmp(current,"long",4)==0)||
             (strncmp(current,"enum",4)==0)||
             (strncmp(current,"goto",4)==0)||
             (strncmp(current,"void",4)==0)||
             (strncmp(current,"else",4)==0)
             )
            ){
          // not valid, skip
          ;
        }else if((currentIdx==5)&&
            (
             (strncmp(current,"break",5)==0)||
             (strncmp(current,"const",5)==0)||
             (strncmp(current,"short",5)==0)||
             (strncmp(current,"float",5)==0)||
             (strncmp(current,"union",5)==0)||
             (strncmp(current,"while",5)==0)
             )
            ){
          // not valid, skip
          ;
        }else if((currentIdx==6)&&
            (
             (strncmp(current,"double",6)==0)||
             (strncmp(current,"return",6)==0)||
             (strncmp(current,"extern",6)==0)||
             (strncmp(current,"signed",6)==0)||
             (strncmp(current,"sizeof",6)==0)||
             (strncmp(current,"static",6)==0)||
             (strncmp(current,"struct",6)==0)||
             (strncmp(current,"switch",6)==0)
             )
            ){
          // not valid, skip
          ;
        }else if((currentIdx==7)&&
            (
             (strncmp(current,"default",7)==0)||
             (strncmp(current,"typedef",7)==0)
             )
            ){
          // not valid, skip
          ;
        }else if((currentIdx==8)&&
            (
             (strncmp(current,"continue",8)==0)||
             (strncmp(current,"register",8)==0)||
             (strncmp(current,"unsigned",8)==0)||
             (strncmp(current,"volatile",8)==0)
             )
            ){
          // not valid, skip
          ;
        }else{
          // check whther every char is valid
          valid=1;
          for(i=0;i<currentIdx;i++){
            if((current[i]=='_')||
                ((current[i]<='z')&&(current[i]>='a'))||
                ((current[i]<='Z')&&(current[i]>='A'))||
                ((current[i]<='9')&&(current[i]>='0'))){
              ;
            }else{
              valid=0;
              break;
            }
          }
          if(valid==1){
            printf("%s\n",current);
            count++;
          }
        }
        inStringFlag=0;
        currentIdx=0;
      }
      if(input[inputIdx]==0){
        printf("PROGRAM OUTPUT:%d identifiers!\n",count);
        exit(0);
      }
    }else{
      inStringFlag=1;
      current[currentIdx]=input[inputIdx];
      currentIdx++;
    }
    inputIdx++;
  }
}

