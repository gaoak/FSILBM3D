#include<unistd.h>
void myfork_(int* p) {
    *p = fork();
}

void mysleep_(int* t) {
    sleep((unsigned int) (*t));
}