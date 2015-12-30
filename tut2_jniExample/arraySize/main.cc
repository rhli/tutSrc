#include <cstdio>
#include <cstdlib>

int calculateSize(int* a){
    printf("%d\n",sizeof(*a)/sizeof(int));
}


int main(){
    int a[7];
    calculateSize(a);
    return 0;
}
