#ifndef HEADERFILEFILEIO_HPP
#define HEADERFILEFILEIO_HPP
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>
static float EOHMARKER = 357.0f;

static float ZONEMARKER = 299.0f;

template <typename T> size_t FindMax(size_t N, const T *data) {
  size_t ind = 0;
  T tmp = data[0];
  for (size_t i = 1; i < N; ++i) {
    if (tmp < data[i]) {
      ind = i;
      tmp = data[i];
    }
  }
  return ind;
}

template <typename T> size_t FindMin(size_t N, const T *data) {
  size_t ind = 0;
  T tmp = data[0];
  for (size_t i = 1; i < N; ++i) {
    if (tmp > data[i]) {
      ind = i;
      tmp = data[i];
    }
  }
  return ind;
}

int BinaryWrite(std::ofstream &ofile, std::string str) {
  int tmp = 0;
  for (size_t i = 0; i < str.size(); ++i) {
    tmp = str[i];
    ofile.write((char *)&tmp, 4);
  }
  tmp = 0;
  ofile.write((char *)&tmp, 4);
  return 4;
}

template <typename T>
int OutputTec360_binary(const std::string filename,
                        const std::vector<std::string> &variables,
                        const std::vector<size_t> &rawN,
                        const std::vector<std::vector<T>> &data) {
  int isdouble = sizeof(T) / 8;
  std::vector<size_t> N = rawN;
  for (size_t i = N.size(); i < 3; ++i) {
    N.push_back(1);
  }
  std::ofstream odata;
  odata.open(filename, std::ios::binary);
  if (!odata.is_open()) {
    printf("error unable to open file %s\n", filename.c_str());
    return -1;
  }
  char tecplotversion[] = "#!TDV112";
  odata.write((char *)tecplotversion, 8);
  int value1 = 1;
  odata.write((char *)&value1, 4);
  int filetype = 0;
  odata.write((char *)&filetype, 4);
  // read file title and variable names
  std::string filetitle = "";
  BinaryWrite(odata, filetitle);
  int nvar = variables.size();
  odata.write((char *)&nvar, 4); // number of variables
  std::vector<std::string> vartitle;
  for (int i = 0; i < nvar; ++i) {
    BinaryWrite(odata, variables[i]);
  }
  float zonemarker299I = ZONEMARKER;
  odata.write((char *)&zonemarker299I, 4);
  // zone title
  std::string zonetitle("ZONE 0");
  BinaryWrite(odata, zonetitle);
  int parentzone = -1;
  odata.write((char *)&parentzone, 4);
  int strandid = -1;
  odata.write((char *)&strandid, 4);
  double soltime = 0.0;
  odata.write((char *)&soltime, 8);
  int unused = -1;
  odata.write((char *)&unused, 4);
  int zonetype = 0;
  odata.write((char *)&zonetype, 4);
  int zero = 0;
  odata.write((char *)&zero, 4);
  odata.write((char *)&zero, 4);
  odata.write((char *)&zero, 4);
  for (int i = 0; i < 3; ++i) {
    int tmp = N[i];
    odata.write((char *)&tmp, 4);
  }

  odata.write((char *)&zero, 4);
  float eohmarker357 = EOHMARKER;
  odata.write((char *)&eohmarker357, 4);
  float zonemarker299II = ZONEMARKER;
  odata.write((char *)&zonemarker299II, 4);
  std::vector<int> binarydatatype(nvar, 1 + (isdouble > 0));
  odata.write((char *)binarydatatype.data(), 4 * nvar);
  odata.write((char *)&zero, 4);
  odata.write((char *)&zero, 4);
  int minus1 = -1;
  odata.write((char *)&minus1, 4);

  size_t datanumber, datasize;
  datanumber = N[0] * N[1] * N[2];
  datasize = N[0] * N[1] * N[2] * sizeof(T);
  for (int i = 0; i < nvar; ++i) {
    double minv = 0., maxv = 1.;
    maxv = data[i][FindMax<T>(datanumber, data[i].data())];
    minv = data[i][FindMin<T>(datanumber, data[i].data())];
    odata.write((char *)&minv, 8);
    odata.write((char *)&maxv, 8);
  }

  for (int i = 0; i < nvar; ++i) {
    odata.write((char *)data[i].data(), datasize);
  }
  odata.close();
  return 0;
}
#endif
