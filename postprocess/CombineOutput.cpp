#include "FileIO.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

typedef float INREAL;
typedef float OUTREAL;
void readGridFiles(const std::string &igridfile, std::vector<OUTREAL> &cData,
                   std::vector<OUTREAL> &dcData);

namespace fs = std::filesystem;

void SubsExactSolution(const std::vector<OUTREAL> &x,
                       const std::vector<OUTREAL> &y,
                       const std::vector<OUTREAL> &z, std::vector<OUTREAL> &u,
                       std::vector<OUTREAL> &v, std::vector<OUTREAL> &w,
                       std::vector<OUTREAL> &p) {
  for (size_t i = 0; i < z.size(); ++i) {
    u[i] -= exp(-0.7071067811865*z[i])*cos(1*12000*0.00628319-0.7071067811865*z[i]);
    v[i] -= 0.;
    w[i] -= 0.;
    p[i] -= 0.;
  }
}

void Integrate(const std::vector<OUTREAL> &dxc, const std::vector<OUTREAL> &dyc,
               const std::vector<OUTREAL> &dzc, const std::vector<OUTREAL> &var,
               int n, OUTREAL &resn, OUTREAL &resinfty) {
  int Nx = dxc.size(), Ny = dyc.size(), Nz = dzc.size();
  int count = 0, cnt2 = 0;
  resinfty = -1;
  // integrate along x
  count = 0;
  cnt2 = 0;
  std::vector<OUTREAL> res1(Ny * Nz, 0.);
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i, ++cnt2) {
        OUTREAL value = fabs(var[cnt2]);
        res1[count] += pow(value, n) * dxc[i];
        if (resinfty < value) {
          resinfty = value;
        }
      }
      ++count;
    }
  }
  // integrate along y
  count = 0;
  cnt2 = 0;
  std::vector<OUTREAL> res2(Nz, 0.);
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j, ++cnt2) {
      res2[count] += res1[cnt2] * dyc[j];
    }
    ++count;
  }
  // integrate along z
  OUTREAL result = 0.;
  for (int k = 0; k < Nz; ++k) {
    result += res2[k] * dzc[k];
  }
  resn = pow(result, 1. / n);
}

void Integrate(const std::vector<OUTREAL> &dxc, const std::vector<OUTREAL> &dyc,
               const std::vector<OUTREAL> &dzc,
               const std::vector<std::vector<OUTREAL>> &vars, int n,
               std::vector<OUTREAL> &resn, std::vector<OUTREAL> &resinfty) {
  resn.resize(vars.size(), 0.);
  resinfty.resize(vars.size(), 0.);
  std::vector<std::vector<OUTREAL>> data = vars;
  SubsExactSolution(data[0], data[1], data[2], data[3], data[4], data[5],
                    data[6]);
  for (size_t v = 0; v < data.size(); ++v) {
    Integrate(dxc, dyc, dzc, data[v], n, resn[v], resinfty[v]);
  }
}

void printErrors(const std::string xgrid_file, const std::string ygrid_file,
                 const std::string zgrid_file,
                 const std::vector<std::string> &vars,
                 const std::vector<std::vector<OUTREAL>> &Stacks) {
  std::vector<OUTREAL> dxc;
  std::vector<OUTREAL> dyc;
  std::vector<OUTREAL> dzc;
  std::vector<OUTREAL> coord;

  readGridFiles(xgrid_file, coord, dxc);
  readGridFiles(ygrid_file, coord, dyc);
  readGridFiles(zgrid_file, coord, dzc);

  std::vector<OUTREAL> resn;
  std::vector<OUTREAL> resinfty;
  Integrate(dxc, dyc, dzc, Stacks, 2, resn, resinfty);
  printf("============================\n");
  for (size_t v = 0; v < Stacks.size(); ++v) {
    printf("VAR = %s, L2 norm %24.16e\n", vars[v].c_str(), resn[v]);
  }
  printf("============================\n");
  for (size_t v = 0; v < Stacks.size(); ++v) {
    printf("VAR = %s, LInfinity norm %24.16e\n", vars[v].c_str(), resinfty[v]);
  }
}

// Used to compare file names containing process numbers
bool compare_filenames_by_process_number(const std::filesystem::path &a,
                                         const std::filesystem::path &b) {
  std::string filename_a = a.filename().string();
  std::string filename_b = b.filename().string();

  // Assuming the file name format is "Flow_00<process number>. dat"
  size_t dot_pos_a = filename_a.find('.');
  size_t dot_pos_b = filename_b.find('.');

  if (dot_pos_a == std::string::npos || dot_pos_b == std::string::npos) {
    // If '.' cannot be found, Then it is considered that the file does not
    // conform to the format and placed at the end of the sorting
    return false;
  }

  std::string process_number_a = filename_a.substr(6, dot_pos_a - 6);
  std::string process_number_b = filename_b.substr(6, dot_pos_b - 6);

  // Convert strings to integers for comparison
  int num_a = std::stoi(process_number_a);
  int num_b = std::stoi(process_number_b);

  return num_a < num_b;
}

void readGridFiles(const std::string &igridfile, std::vector<OUTREAL> &cData,
                   std::vector<OUTREAL> &dcData) {
  std::ifstream gridfile(igridfile);
  std::string line;
  std::vector<OUTREAL> Data;
  cData.clear();
  dcData.clear();
  if (!gridfile.is_open()) {
    std::cerr << "Failed to open file: " << igridfile << std::endl;
    return;
  }

  int id;
  OUTREAL x;
  while (std::getline(
      gridfile,
      line)) { // If getline successfully reads a line, it will return true
    std::istringstream iss(line);
    iss >> id >> x;
    if (iss.good()) {
      Data.push_back(x);
    }
  }
  for (int i = 0; i < Data.size() - 1; ++i) {
    cData.push_back(0.5 * (Data[i] + Data[i + 1]));
    dcData.push_back(Data[i + 1] - Data[i]);
  }
  gridfile.close();
}

void readGridFiles(const std::string &igridfile, std::vector<OUTREAL> &cData) {
  std::vector<OUTREAL> dcData;
  readGridFiles(igridfile, cData, dcData);
}

void ExtractIb(std::vector<char> ibData, OUTREAL *data, int number) {
  for (int i = 0; i < number; ++i) {
    int j = i / 4;
    int k = (i % 4) * 2;
    data[i] = int((ibData[j] >> k) & 3);
  }
}

void stackDataFromFiles(const std::vector<fs::path> &files,
                        const std::string xgrid_file,
                        const std::string ygrid_file,
                        const std::string zgrid_file,
                        std::vector<std::vector<OUTREAL>> &Stacks, int &NXc,
                        int &NYc, int &NZc, bool showIb) {
  std::vector<OUTREAL> xcData;
  std::vector<OUTREAL> ycData;
  std::vector<OUTREAL> zcData;
  int zc_start, zc_end;

  readGridFiles(xgrid_file, xcData);
  readGridFiles(ygrid_file, ycData);
  readGridFiles(zgrid_file, zcData);
  NXc = xcData.size();
  NYc = ycData.size();
  NZc = zcData.size();
  if (showIb) {
    Stacks.resize(8);
  } else {
    Stacks.resize(7);
  }
  for (size_t i = 0; i < Stacks.size(); ++i) {
    Stacks[i].resize(NXc * NYc * NZc);
  }

  int count = 0;
  for (int k = 0; k < NZc; ++k) {
    for (int j = 0; j < NYc; ++j) {
      for (int i = 0; i < NXc; ++i) {
        Stacks[0][count] = xcData[i];
        Stacks[1][count] = ycData[j];
        Stacks[2][count] = zcData[k];
        ++count;
      }
    }
  }
  // Read files in sorted order
  std::set<int> zids;
  for (int i = 1; i <= NZc; ++i) {
    zids.insert(i);
  }
  for (const auto &file : files) {
    std::string filename = file.string();
    std::ifstream infile(filename, std::ios::binary);
    if (!infile.is_open()) {
      std::cerr << "Failed to open file: " << filename << std::endl;
      return;
    }
    infile.read(reinterpret_cast<char *>(&zc_start), sizeof(zc_start));
    infile.read(reinterpret_cast<char *>(&zc_end), sizeof(zc_end));
    // Read data
    int offset = NXc * NYc * (zc_start - 1);
    int ndata = NXc * NYc * (zc_end - zc_start + 1);
    int datasize = sizeof(INREAL) * ndata;
    if (std::is_same<INREAL, OUTREAL>::value) {
      infile.read((char *)&Stacks[3][offset], datasize);
      infile.read((char *)&Stacks[4][offset], datasize);
      infile.read((char *)&Stacks[5][offset], datasize);
      infile.read((char *)&Stacks[6][offset], datasize);
    } else {
      std::vector<INREAL> tempBuffer(datasize / sizeof(INREAL));
      for (int i = 0; i < 4; ++i) {
        infile.read(reinterpret_cast<char *>(tempBuffer.data()), datasize);
        for (size_t j = 0; j < tempBuffer.size(); j++) {
          Stacks[i + 3][offset + j] = static_cast<OUTREAL>(tempBuffer[j]);
        }
      }
    }
    if (showIb) {
      datasize = std::ceil(ndata / 4);
      std::vector<char> ibdata(datasize);
      infile.read((char *)&ibdata[0], datasize);
      ExtractIb(ibdata, &Stacks[7][offset], ndata);
    }
    infile.close();
    for (int i = zc_start; i <= zc_end; ++i)
      zids.erase(i);
  }
  if (!zids.empty()) {
    std::cout << "Inconsistent data size " << std::endl;
    std::cout << "missling z-slices: ";
    for (auto i : zids) {
      std::cout << i << " ";
    }
    std::cout << "." << std::endl;
  }
}

int main(int argc, char *argv[]) {
  if (argc < 6) {
    std::cout << "Error, insufficient input parameters, Usage: Combine "
                 "xgrid.dat ygrid.dat zgrid.dat Flow/Time_0000 0.plt"
              << std::endl;
    exit(-1);
  }
  std::string optIblankStr("--iblank");
  bool showIb = false;
  for (int i = 1; i < argc; ++i) {
    if (optIblankStr.compare(argv[i]) == 0) {
      showIb = true;
      break;
    }
  }
  std::string optErrorStr("--error");
  bool showError = false;
  for (int i = 1; i < argc; ++i) {
    if (optErrorStr.compare(argv[i]) == 0) {
      showError = true;
      break;
    }
  }
  std::string xgrid_file(argv[1]);
  std::string ygrid_file(argv[2]);
  std::string zgrid_file(argv[3]);
  fs::path entry_path(argv[4]);
  std::string plt_file(argv[5]);
  std::vector<int> rawN(3);
  std::vector<std::vector<OUTREAL>> Stacks;
  std::vector<std::string> var = {"x", "y", "z", "u", "v", "w", "p"};
  if (showIb) {
    var.push_back("Ib");
  }
  std::vector<fs::path> files;
  // Collect all eligible file paths
  for (const auto &file : fs::directory_iterator(entry_path)) {
    if (file.is_regular_file() &&
        file.path().filename().string().find("Flow_") == 0) {
      files.push_back(file.path());
    }
  }
  // Sort file paths based on process numbers
  std::sort(files.begin(), files.end(), compare_filenames_by_process_number);
  stackDataFromFiles(files, xgrid_file, ygrid_file, zgrid_file, Stacks, rawN[0],
                     rawN[1], rawN[2], showIb);
  OutputTec360_binary(plt_file, var, rawN, Stacks);
  std::cout << "Write " << plt_file << std::endl;
  if (showError) {
    printErrors(xgrid_file, ygrid_file, zgrid_file, var, Stacks);
  }
  return 0;
}
