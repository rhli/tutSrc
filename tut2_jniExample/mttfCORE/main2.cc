# include <stdio.h>
# include <stdlib.h>

# define N 10
# define K 8

int main(int argc, char** argv){
    int stateNum=N-K+2;
    printf("1/(");
    int indexer =0 ;
    for(int i=N-K+1;i>=0;i--){
        if(i==N-K+1){
            printf("1");
        }else if(i==N-K){
            printf("+s/(%d*l)",N-K);
        }else{
            printf("+");
            for(int j=0;j<indexer-1;j++){
                printf("(s+%d*l+u)*",N-K+j);
            }
            printf("s/(");
            for(int j=0;j<indexer;j++){
                if(j!=indexer-1){
                    printf("%d*l*",N-K+j);
                }else{
                    printf("%d*l",N-K+j);
                }
            }
            printf(")");
        }
        indexer ++;
    }
    printf(")\n");

    /*
     * CORE
     */
    printf("CORE\n");
    printf("1/(");
    indexer =0 ;
    for(int i=N-K+1;i>=0;i--){
        if(i==N-K+1){
            printf("1");
        }else if(i==N-K){
            printf("+s/(%d*l)",N-K);
        }else{
            printf("+");
            for(int j=0;j<indexer-1;j++){
                printf("(s+%d*l+%lf*bsr)*",N-K+j,(double)3600*24*365/1000/8/K);
            }
            printf("s/(");
            for(int j=0;j<indexer;j++){
                if(j!=indexer-1){
                    printf("%d*l*",N-K+j);
                }else{
                    printf("%d*l",N-K+j);
                }
            }
            printf(")");
        }
        indexer ++;
    }
    printf(")\n");
    return 0;
}
