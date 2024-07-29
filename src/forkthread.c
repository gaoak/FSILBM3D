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