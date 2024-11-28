#include<cstdint>
#include<map>
#include<vector>
#include<iostream>
extern "C" {
    void initumap(int32_t &np);
    void thread_setumap(int32_t &p, int64_t &index, double &value);
    void mergeumap();
    void setumap(int64_t &index, double &value);
    void getumap(int64_t &index, double &value);
    void printumap();
}

std::vector<std::map<int64_t, double>> threadmap;
std::map<int64_t, double> umap;

void initumap(int32_t &np) {
    threadmap.resize(np);
    umap.clear();
}

void thread_setumap(int32_t &p, int64_t &index, double &value) {
    threadmap[p][index] = value;
}

void mergeumap() {
    for(auto &m : threadmap) {
        for(const auto &it : m) {
            umap[it.first] = it.second;
        }
        m.clear();
    }
    threadmap.clear();
}

void setumap(int64_t &index, double &value) {
    umap[index] = value;
}

void getumap(int64_t &index, double &value) {
    value = umap[index];
}


void printumap() {
    std::cout << "umap with size " << umap.size() << std::endl;
    for(const auto &it:umap) {
        std::cout << it.first << ", " << it.second << std::endl;
    }
    std::cout << "================" << std::endl;
}