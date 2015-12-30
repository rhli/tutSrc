# include <stdio.h>
# include <stdlib.h>
# include <string.h>

int main(int argc,char** argv){
    int index=atoi(argv[3]);
    printf("(");
    printf("(");
    printf("(s+u+%d*l)*%s-u*%s",19-index,argv[2],argv[1]);
    printf(")");
    printf("/");
    printf("(");
    printf("%d*l",20-index);
    printf(")");
    printf(")");
    printf("\n");
    printf("\n");
    printf("\n");
    printf("\\(");
    printf("\\(");
    printf("\\(s+u+%d*l\\)*",19-index);
    for(int i=0;i<strlen(argv[2]);i++){
        if((argv[2][i]=='(')||(argv[2][i]==')')){
            printf("\\%c",argv[2][i]);
        }else{
            printf("%c",argv[2][i]);
        }
    }
    printf("-u*");
    for(int i=0;i<strlen(argv[1]);i++){
        if((argv[1][i]=='(')||(argv[1][i]==')')){
            printf("\\%c",argv[1][i]);
        }else{
            printf("%c",argv[1][i]);
        }
    }
    printf("\\)");
    printf("/");
    printf("\\(");
    printf("%d*l",20-index);
    printf("\\)");
    printf("\\)");
    printf("\n");
    return 0;
}
