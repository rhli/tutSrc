# include <stdio.h>
# include <stdlib.h>

//# define LAMBDA 0.06
//# define B 0.1 //Gbps
//# define S 4 //TB

int main(int argc, char** argv){
    if(argc!=7){
        printf("Usage: %s (n) (k) (relayer or OBO) (core or RS) (bsrate) (lambda)\n",argv[0]);
        exit(0);
    }
    int n=atoi(argv[1]);
    int k=atoi(argv[2]);
    /*
     * 1 for relayer
     * 2 for one-by-one
     */
    int type=atoi(argv[3]);
    /*
     * 1 for CORE
     * 2 for erasure code
     */
    int code=atoi(argv[4]);
    double bsrate = atof(argv[5]);
    double LAMBDA = atof(argv[6]);
    int stateNum=n-k+2;
    double* transitionMat=(double*)calloc(stateNum*stateNum,sizeof(double));
    if(type==1){
        if(code==1){
            for(int i=0;i<stateNum;i++){
                if(i!=0 && i!=stateNum-1){ 
                    transitionMat[i*stateNum]
                        =(double)bsrate*(n-k)/i/(n-i)*3600*24*365/1000/8;
                }
                if(i!=stateNum-1){ 
                    transitionMat[i*stateNum+i+1]
                        =((double)(n-i))*LAMBDA;
                }
            }
        }else if(code==2){
            for(int i=0;i<stateNum;i++){
                if(i!=0 && i!=stateNum-1){ 
                    transitionMat[i*stateNum]
                        =(double)bsrate/k*3600*24*365/1000/8;
                }
                if(i!=stateNum-1){ 
                    transitionMat[i*stateNum+i+1]
                        =((double)(n-i))*LAMBDA;
                }
            }
        }
        for(int i=0;i<stateNum;i++){
            for(int j=0;j<stateNum;j++){
                if(j!=i){
                    transitionMat[i*stateNum+i]-=transitionMat[i*stateNum+j];
                }
            }
        }
    }else if(type==2){
        if(code==1){
            for(int i=0;i<stateNum;i++){
                if(i!=0 && i!=stateNum-1){ 
                    transitionMat[i*stateNum+i-1]
                        =(double)bsrate/k*3600*24*365/1000/8;
                }
                if(i!=stateNum-1){ 
                    transitionMat[i*stateNum+i+1]
                        =((double)(n-i))*LAMBDA;
                }
            }
            for(int i=0;i<stateNum;i++){
                for(int j=0;j<stateNum;j++){
                    if(j!=i){
                        transitionMat[i*stateNum+i]-=transitionMat[i*stateNum+j];
                    }
                }
            }
        }else if(code==2){
            for(int i=0;i<stateNum;i++){
                if(i!=0 && i!=stateNum-1){ 
                    transitionMat[i*stateNum+i-1]
                        =(double)bsrate/k*3600*24*365/1000/8;
                }
                if(i!=stateNum-1){ 
                    transitionMat[i*stateNum+i+1]
                        =((double)(n-i))*LAMBDA;
                }
            }
            for(int i=0;i<stateNum;i++){
                for(int j=0;j<stateNum;j++){
                    if(j!=i){
                        transitionMat[i*stateNum+i]-=transitionMat[i*stateNum+j];
                    }
                }
            }
        }
    }
    printf("Transition Matrix:\n");
    for(int i=0;i<stateNum;i++){
        for(int j=0;j<stateNum;j++){
            printf("%10lf ",transitionMat[i*stateNum+j]);
        }
        printf("\n");
    }
    for(int i=0;i<stateNum-1;i++){
        printf("%lf,",transitionMat[i*stateNum]);
    }
    printf("\n");
    printf("[");
    for(int i=0;i<stateNum;i++){
        printf("x%d",i);
        if(i!=stateNum-1){
            printf(",");
        }
    }
    printf("]=");
    printf("solve(");
    for(int i=0;i<stateNum;i++){
        if(i==0){
            printf("'s*x0-1=");
        }else{
            printf("'s*x%d=",i);
        }
        int index=0;
        for(int j=0;j<stateNum;j++){
            if(index==0){
                if(transitionMat[j*stateNum+i]==0) continue;
                printf("(%lf)*x%d",transitionMat[j*stateNum+i],j);
                index++;
            }else{
                if(transitionMat[j*stateNum+i]==0) continue;
                printf("+(%lf)*x%d",transitionMat[j*stateNum+i],j);
                index++;
            }
        }
        if(i==stateNum-1){
            printf("',");
        }else{
            printf("',");
        }
    }
    for(int i=0;i<stateNum;i++){
        if(i==stateNum-1){
            printf("'x%d'",i);
        }else{
            printf("'x%d',",i);
        }
    }
    printf(")");
    //printf("}]\n");
    printf("\n");
    printf("[");
    for(int i=0;i<stateNum;i++){
        printf("x%d",i);
        if(i!=stateNum-1){
            printf(",");
        }
    }
    printf("]=solve(");
    if(type==1){
        if(code==1){
            for(int i=0;i<stateNum;i++){
                if(i==0){
                    printf("'s*x0-1=");
                    for(int j=0;j<stateNum-1;j++){
                        if(j==0){
                            printf("(-%d*l)*x0",n);
                        }else{
                            printf("+u%d*x%d",j,j);
                            //printf("+u*%lf*x%d",(double)(n-k)*k/(n-j)/j,j);
                        }
                    }
                }else if(i<stateNum-1){
                    printf("'s*x%d=%d*l*x%d-(%d*l+u%d)*x%d",i,n-i+1,i-1,n-i,i,i);
                    //printf("'s*x%d=%d*l*x%d-(%d*l+u*%lf)*x%d",
                    //        i,n-i+1,i-1,n-i,(double)(n-k)*k/(n-i)/i,i);
                }else{
                    printf("'s*x%d=%d*l*x%d",i,n-i+1,i-1);
                }
                printf("',");
            }
        }else if(code==2){
            for(int i=0;i<stateNum;i++){
                if(i==0){
                    printf("'s*x0-1=");
                    for(int j=0;j<stateNum-1;j++){
                        if(j==0){
                            printf("(-%d*l)*x0",n);
                        }else{
                            printf("+u*x%d",j);
                        }
                    }
                }else if(i<stateNum-1){
                    printf("'s*x%d=%d*l*x%d-(%d*l+u)*x%d",i,n-i+1,i-1,n-i,i);
                }else{
                    printf("'s*x%d=%d*l*x%d",i,n-i+1,i-1);
                }
                printf("',");
            }
        }
    }else{
        for(int i=0;i<stateNum;i++){
            if(i==0){
                printf("'s*x0-1=-%d*l*x0+u1*x1',",n);
            }else if(i<stateNum-2){
                printf("'s*x%d=%d*l*x%d-(u%d+%d*l)*x%d+u%d*x%d',",
                        i,n-i+1,i-1,i,n-i,i,i+1,i+1);
            }else if(i==stateNum-2){
                printf("'s*x%d=%d*l*x%d-(u%d+%d*l)*x%d',",i,n-i+1,i-1,i,n-i,i);
            }else{
                printf("'s*x%d=%d*l*x%d',",i,i-2,i-1);
            }
        }
    }
    for(int i=0;i<stateNum;i++){
        if(i==stateNum-1){
            printf("\'x%d\')",i);
        }else{
            printf("\'x%d\',",i);
        }
    }
    printf("\n");
    return 0;
}




