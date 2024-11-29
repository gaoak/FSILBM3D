#include<cstdint>
#include<map>
#include<vector>
#include<iostream>
extern "C" {
    void initumap(int32_t &np);
    void thread_setumap(int32_t &p, int64_t &index, std::vector<double> &value);
    void mergeumap();
    void setumap(int64_t &index, std::vector<double> &value);
    void getumap(int64_t &index, std::vector<double> &value);
    void printumap();
}

std::vector<std::map<int64_t, std::vector<double>>> threadmap;
std::map<int64_t, std::vector<double>> umap;

void initumap(int32_t &np) {
    threadmap.resize(np);
    umap.clear();
    for(auto &m : threadmap) {
        m.clear();
    }
}

void thread_setumap(int32_t &p, int64_t &index, std::vector<double> &value) {
    threadmap[p][index] = value;
}

void mergeumap() {
    for(auto &m : threadmap) {
        for(const auto &it : m) {
            umap[it.first] = it.second;
        }
        m.clear();
    }
}

void setumap(int64_t &index, std::vector<double> &value) {
    umap[index] = value;
}

void getumap(int64_t &index, std::vector<double> &value) {
    value = umap[index];
}


void printumap() {
    std::cout << "umap with size " << umap.size() << std::endl;
    for(const auto &it:umap) {
        std::cout << it.first << ", [";
        for (size_t i = 0; i < it.second.size(); ++i) {
            std::cout << it.second[i];
            if (i < it.second.size() - 1) {
                std::cout << ", ";  // 添加逗号分隔符
            }
        }
        std::cout << "]" << std::endl;
    }
    std::cout << "================" << std::endl;
}