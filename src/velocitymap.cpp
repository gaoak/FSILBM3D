#include<cstdint>
#include<set>
#include<map>
#include<vector>
#include<iostream>
extern "C" {
    void initumap(int32_t &np);
    void thread_adduindex(int32_t &p, int64_t &index);
    void mergeumap();
    void allocateuarray();
    void setumap(int64_t &index, double* value);
    void getumap(int64_t &index, double* value);
    void printumap();
}

std::vector<std::set<int64_t>> threadset;
std::map<int64_t, int> umap;
std::vector<double> velocity;

void initumap(int32_t &np) {
    threadset.resize(np);
    umap.clear();
    for(auto &m : threadset) {
        m.clear();
    }
}

void thread_adduindex(int32_t &p, int64_t &index) {
    threadset[p].insert(index);
}

void mergeumap() {
    int count = int(umap.size());
    for(auto &m : threadset) {
        for(const auto &it : m) {
            umap[it] = 3*count;
            ++count;
        }
        m.clear();
    }
}

void allocateuarray() {
    velocity.resize(3*umap.size());
}

void setumap(int64_t &index, double* value) {
    int i = umap[index];
    velocity[i  ] = value[0];
    velocity[i+1] = value[1];
    velocity[i+2] = value[2];
}

void getumap(int64_t &index, double value[3]) {
    int i = umap[index];
    value[0] = velocity[i  ];
    value[1] = velocity[i+1];
    value[2] = velocity[i+2];
}

void printumap() {
    std::cout << "umap with size " << umap.size() << std::endl;
    for(const auto &it:umap) {
        std::cout << it.first << ", [";
        for (size_t i = 0; i < 3; ++i) {
            std::cout << velocity[it.second+i];
            if (i < 2) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
    }
    std::cout << "================" << std::endl;
}