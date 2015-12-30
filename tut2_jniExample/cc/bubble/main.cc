#include <stdio.h>
#define N 5

int main(){
  int a[N];
  int *q=NULL;
  void BubbleSort(int *a,int n);
  printf("Please input the numbers:\n");
  for(q=a;q<&a[N];q++)scanf("%d",q);
  BubbleSort(a,N);
  printf("The sorted numbers:\n");
  for(q=a;q<a+N;q++)printf("%d ",*q);
  printf("\n");
}

void BubbleSort(int *a,int n){
  int temp;
  int i;
  int *p=NULL;
  for(i=0;i<n-1;i++){ 
    for(p=a;p<a+n-i-1;p++) 
      if(*p<*(p+1)){ 
        temp=*p;
        *p=*(p+1);
        *(p+1)=temp; 
      }
  }
}

