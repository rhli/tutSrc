#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void reverseString(char*s,int start,int end){
  int i; 
  int pair=(end-start+1)/2; 
  for (i=0;i<pair;i++){ 
    s[start+i]^=s[end-i]; 
    s[end-i]^=s[start+i]; 
    s[start+i]^=s[end-i];
  }
}

void reverseWords(char *s) { 
  int i; 
  int totalLen=strlen(s); 
  int startLoc=0; 
  bool inWord=false; 
  char workspace[totalLen+1]; 
  for(i=0;i<totalLen+1;i++){ 
    workspace[i]=0; 
  } 
  int workspacepos=0; 
  for(i=0;i<totalLen+1;i++){ 
    if((s[i]==' ')||(s[i]==0)){ 
      if(inWord){ 
        strncpy(workspace+workspacepos,s+startLoc,i-startLoc+1); 
        reverseString(workspace,workspacepos,workspacepos+i-startLoc-1); 
        workspacepos=workspacepos+i-startLoc+1;
      } 
      inWord=false; 
    }else{ 
      if(!inWord){ 
        inWord=true; 
        startLoc=i; 
      } 
    } 
  } 
  workspacepos-=2;
  if(workspacepos<0)workspacepos=0;
  reverseString(workspace,0,workspacepos);
  strncpy(s,workspace,workspacepos+2);
}

int main(){
  char s[100];
  gets(s);
  reverseWords(s);
  printf("%s*\n",s);
}


