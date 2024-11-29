#include<cstdint>
#include<set>
#include<map>
#include<vector>
#include<iostream>
#include<iterator>
extern "C" {
    void initumap(int32_t &np);
    void thread_adduindex(int32_t &p, int64_t &index);
    void mergeumap();
    void allocateuarray();
    void setumap(int64_t &index, double* value);
    void addumap(int64_t &index, double* value);
    void getumap(int64_t &index, double* value);
    void printumap();
    void findindex(uint16_t* index, int32_t &i, int32_t &j, int32_t &k);
    void findijk(uint16_t* index, int32_t &i, int32_t &j, int32_t &k);
    void inititerator(int32_t &np, int32_t *ndata, int64_t &index);
    void getiterator(int32_t &p, int64_t &index);
    void nextiterator(int32_t &p);
}

std::vector<std::set<int64_t>> threadset;
std::map<int64_t, int> umap;
std::vector<double> velocity;
std::vector<std::map<int64_t, int>::iterator> iter;
void initumap(int32_t &np) {
    threadset.resize(np);
    umap.clear();
    for(auto &m : threadset) {
        m.clear();
    }
}

void inititerator(int32_t &np, int32_t *ndata, int64_t &index) {
    int32_t ntmp = umap.size() / np;
    for(int p=0; p<np; ++p) {
        ndata[p] = ntmp;
        if(p==np-1) {
            ndata[p] += umap.size() - np * ntmp;
        }
    }
    iter.resize(np);
    iter[0] = umap.begin();
    for(int p=1; p<np; ++p) {
        std::advance(iter[p], ndata[p]);
    }
}

void getiterator(int32_t &p, int64_t &index) {
    index = iter[p]->first;
}

void nextiterator(int32_t &p) {
    ++iter[p];
}

void thread_adduindex(int32_t &p, int64_t &index) {
    threadset[p].insert(index);
}

void mergeumap() {
    int count = int(umap.size());
    for(auto &m : threadset) {
        for(const auto &it : m) {
            umap[it] = count;
            count += 3;
        }
        m.clear();
    }
}

void allocateuarray() {
    velocity.resize(3*umap.size(), 0.);
}

void setumap(int64_t &index, double* value) {
    int i = umap[index];
    velocity[i  ] = value[0];
    velocity[i+1] = value[1];
    velocity[i+2] = value[2];
}

void addumap(int64_t &index, double* value) {
    int i = umap[index];
    velocity[i  ] += value[0];
    velocity[i+1] += value[1];
    velocity[i+2] += value[2];
}

void getumap(int64_t &index, double* value) {
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

void findindex(uint16_t* index, int32_t &i, int32_t &j, int32_t &k) {
    index[0] = i;
    index[1] = j;
    index[2] = k;
    index[3] = 0;
}

void findijk(uint16_t* index, int32_t &i, int32_t &j, int32_t &k) {
    i = index[0];
    j = index[1];
    k = index[2];
}
