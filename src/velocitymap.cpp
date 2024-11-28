#include<cstdint>
#include<map>
#include<iostream>
extern "C" {
    void setumap(int64_t & index, double* value);
    void getumap(int64_t & index, double* value);
    void clearumap();
    void printumap();
}


std::map<int64_t, double> umap;
void setumap(int64_t & index, double* value) {
    umap[index] = *value;
}

void getumap(int64_t & index, double* value) {
    umap[index] = *value;
}

void clearumap() {
    umap.clear();
}

void printumap() {
    std::cout << "umap with size " << umap.size() << std::endl;
    for(const auto &it:umap) {
        std::cout << it.first << ", " << it.second << std::endl;
    }
    std::cout << "================" << std::endl;
}