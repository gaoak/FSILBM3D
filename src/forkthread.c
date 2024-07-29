#include<unistd.h>
#include<stdlib.h>
void myfork_(int* p) {
    *p = fork();
}

void mysleep_(int* t) {
    sleep((unsigned int) (*t));
}

void myexit_(int* p) {
    exit(*p);
}