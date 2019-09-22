#include<stdio.h>
#include<math.h>

int main(int argc, char *argv[]){
    int c;
    char *filename;
    FILE *fp;
   
    if (argc == 2){
        fp = fopen(argv[1],"r");
        if (fp == NULL){
        printf("%s file not open!\n",argv[1]);
        return -1;
        }
    }
    else if (argc < 2 ){
        printf("Please input a specqq file\n");
        fscanf("%s",filename);
        fp = fopen(filename,"r");
        if (fp == NULL){
        printf("%s file not open!\n",argv[1]);
        return -1;
        }
    }
        
    while((c =fgetc(fp)) != EOF && c != "\n");
    while(fscanf(fp,))


    fclose(fp);
	




 


}
