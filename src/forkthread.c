#include<unistd.h>
#include<stdlib.h>
#include<sys/wait.h>
void myfork_(int* p) {
    *p = fork();
}

void mywait_() {
    int* status;
    wait(status);
}

void mysleep_(int* t) {
    sleep((unsigned int) (*t));
}

void myexit_(int* p) {
    exit(*p);
}

//#include<stdio.h>
/* int main() {
    int i;
    for(i=0; i<3; ++i) {
        int p = fork();
        if(p==0) {
            printf("child start\n");
            sleep(1+i);
            printf("child end\n");
            exit(0);
        }
    }
    sleep(1);
    for(i=0; i<3; ++i) {
        mywait_();
    }
    printf("end wait\n");

    for(i=0; i<3; ++i) {
        int p = fork();
        if(p==0) {
            printf("child2 start\n");
            sleep(i+2);
            printf("child2 end\n");
            exit(0);
        }
    }
    sleep(1);
    for(i=0; i<3; ++i) {
        mywait_();
    }
    printf("end2 wait\n");
} */