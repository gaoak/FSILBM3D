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

namespace fs = std::filesystem;

// void SubsExactSolution(const std::vector<OUTREAL> &x,
//                        const std::vector<OUTREAL> &y,
//                        const std::vector<OUTREAL> &z, std::vector<OUTREAL>
//                        &u, std::vector<OUTREAL> &v, std::vector<OUTREAL> &w,
//                        std::vector<OUTREAL> &p) {
//   for (size_t i = 0; i < z.size(); ++i) {
//     u[i] -=
//     exp(-0.7071067811865*z[i])*cos(1*12000*0.00628319-0.7071067811865*z[i]);
//     v[i] -= 0.;
//     w[i] -= 0.;
//     p[i] -= 0.;
//   }
// }

// void Integrate(const std::vector<OUTREAL> &dxc, const std::vector<OUTREAL>
// &dyc,
//                const std::vector<OUTREAL> &dzc, const std::vector<OUTREAL>
//                &var, int n, OUTREAL &resn, OUTREAL &resinfty) {
//   int Nx = dxc.size(), Ny = dyc.size(), Nz = dzc.size();
//   int count = 0, cnt2 = 0;
//   resinfty = -1;
//   // integrate along x
//   count = 0;
//   cnt2 = 0;
//   std::vector<OUTREAL> res1(Ny * Nz, 0.);
//   for (int k = 0; k < Nz; ++k) {
//     for (int j = 0; j < Ny; ++j) {
//       for (int i = 0; i < Nx; ++i, ++cnt2) {
//         OUTREAL value = fabs(var[cnt2]);
//         res1[count] += pow(value, n) * dxc[i];
//         if (resinfty < value) {
//           resinfty = value;
//         }
//       }
//       ++count;
//     }
//   }
//   // integrate along y
//   count = 0;
//   cnt2 = 0;
//   std::vector<OUTREAL> res2(Nz, 0.);
//   for (int k = 0; k < Nz; ++k) {
//     for (int j = 0; j < Ny; ++j, ++cnt2) {
//       res2[count] += res1[cnt2] * dyc[j];
//     }
//     ++count;
//   }
//   // integrate along z
//   OUTREAL result = 0.;
//   for (int k = 0; k < Nz; ++k) {
//     result += res2[k] * dzc[k];
//   }
//   resn = pow(result, 1. / n);
// }

// void Integrate(const std::vector<OUTREAL> &dxc, const std::vector<OUTREAL>
// &dyc,
//                const std::vector<OUTREAL> &dzc,
//                const std::vector<std::vector<OUTREAL>> &vars, int n,
//                std::vector<OUTREAL> &resn, std::vector<OUTREAL> &resinfty) {
//   resn.resize(vars.size(), 0.);
//   resinfty.resize(vars.size(), 0.);
//   std::vector<std::vector<OUTREAL>> data = vars;
//   SubsExactSolution(data[0], data[1], data[2], data[3], data[4], data[5],
//                     data[6]);
//   for (size_t v = 0; v < data.size(); ++v) {
//     Integrate(dxc, dyc, dzc, data[v], n, resn[v], resinfty[v]);
//   }
// }

// void printErrors(const std::string mesh_file,
//                  const std::vector<std::string> &vars,
//                  const std::vector<std::vector<OUTREAL>> &Stacks) {
//   std::vector<OUTREAL> dxc;
//   std::vector<OUTREAL> dyc;
//   std::vector<OUTREAL> dzc;
//   std::vector<OUTREAL> coord;

//   readGridFiles(xgrid_file, coord, dxc);
//   readGridFiles(ygrid_file, coord, dyc);
//   readGridFiles(zgrid_file, coord, dzc);

//   std::vector<OUTREAL> resn;
//   std::vector<OUTREAL> resinfty;
//   Integrate(dxc, dyc, dzc, Stacks, 2, resn, resinfty);
//   printf("============================\n");
//   for (size_t v = 0; v < Stacks.size(); ++v) {
//     printf("VAR = %s, L2 norm %24.16e\n", vars[v].c_str(), resn[v]);
//   }
//   printf("============================\n");
//   for (size_t v = 0; v < Stacks.size(); ++v) {
//     printf("VAR = %s, LInfinity norm %24.16e\n", vars[v].c_str(),
//     resinfty[v]);
//   }
// }

// Used to compare file names containing process numbers
// bool compare_filenames_by_process_number(const std::filesystem::path &a,
//                                          const std::filesystem::path &b) {
//   std::string filename_a = a.filename().string();
//   std::string filename_b = b.filename().string();

//   // Assuming the file name format is "Flow_00<process number>. dat"
//   size_t dot_pos_a = filename_a.find('.');
//   size_t dot_pos_b = filename_b.find('.');

//   if (dot_pos_a == std::string::npos || dot_pos_b == std::string::npos) {
//     // If '.' cannot be found, Then it is considered that the file does not
//     // conform to the format and placed at the end of the sorting
//     return false;
//   }

//   std::string process_number_a = filename_a.substr(6, dot_pos_a - 6);
//   std::string process_number_b = filename_b.substr(6, dot_pos_b - 6);

//   // Convert strings to integers for comparison
//   int num_a = std::stoi(process_number_a);
//   int num_b = std::stoi(process_number_b);

//   return num_a < num_b;
// }

void readGridFiles(const std::string &igridfile, std::vector<OUTREAL> &xcData,
                   std::vector<OUTREAL> &ycData, std::vector<OUTREAL> &zcData) {
  std::ifstream gridfile(igridfile);
  std::string line;
  std::vector<OUTREAL> Data;
  xcData.clear();
  ycData.clear();
  zcData.clear();
  if (!gridfile.is_open()) {
    std::cerr << "Failed to open file: " << igridfile << std::endl;
    return;
  }

  int id;
  int num = 0;
  int gridNum, dimension = 1;
  OUTREAL x;
  while (std::getline(
      gridfile,
      line)) { // If getline successfully reads a line, it will return true
    std::istringstream iss(line);
    if (num == 0) {
      iss >> gridNum;
    } else {
      iss >> id >> x;
      if (iss.good()) {
        Data.push_back(x);
      }
    }
    num++;
    if (num == gridNum + 1) {
      if (dimension == 1) {
        xcData = Data;
      } else if (dimension == 2) {
        ycData = Data;
      } else {
        zcData = Data;
      }
      Data.clear();
      dimension++;
      num = 0;
    }
  }
  gridfile.close();
}

// void ExtractIb(std::vector<char> ibData, OUTREAL *data, int number) {
//   for (int i = 0; i < number; ++i) {
//     int j = i / 4;
//     int k = (i % 4) * 2;
//     data[i] = int((ibData[j] >> k) & 3);
//   }
// }

void stackDataFromFiles(const std::string filename,
                        std::vector<std::vector<OUTREAL>> &Stacks, size_t &NXc,
                        size_t &NYc, size_t &NZc) {
  std::vector<std::vector<OUTREAL>> tempStacks;
  std::ifstream infile(filename, std::ios::binary);
  if (!infile.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }
  int nx, ny, nz, id;
  infile.read(reinterpret_cast<char *>(&nx), sizeof(nx));
  infile.read(reinterpret_cast<char *>(&ny), sizeof(ny));
  infile.read(reinterpret_cast<char *>(&nz), sizeof(nz));
  infile.read(reinterpret_cast<char *>(&id), sizeof(id));
  double xmin, ymin, zmin, dh;
  infile.read(reinterpret_cast<char *>(&xmin), sizeof(xmin));
  infile.read(reinterpret_cast<char *>(&ymin), sizeof(ymin));
  infile.read(reinterpret_cast<char *>(&zmin), sizeof(zmin));
  infile.read(reinterpret_cast<char *>(&dh), sizeof(dh));

  NXc = nx;
  NYc = ny;
  NZc = nz;
  Stacks.resize(6);
  for (size_t i = 0; i < Stacks.size(); ++i) {
    Stacks[i].resize(NXc * NYc * NZc);
  }
  tempStacks.resize(3);
  for (size_t i = 0; i < tempStacks.size(); ++i) {
    tempStacks[i].resize(NXc * NYc * NZc);
  }

  // Read data
  size_t offset = 0;
  size_t ndata = NXc * NYc * NZc;
  size_t datasize = sizeof(INREAL) * ndata;
  if (std::is_same<INREAL, OUTREAL>::value) {
    infile.read((char *)&tempStacks[0][offset], datasize);
    infile.read((char *)&tempStacks[1][offset], datasize);
    infile.read((char *)&tempStacks[2][offset], datasize);
  } else {
    std::vector<INREAL> tempBuffer(datasize / sizeof(INREAL));
    for (size_t i = 0; i < 3; ++i) {
      infile.read(reinterpret_cast<char *>(tempBuffer.data()), datasize);
      for (size_t j = 0; j < tempBuffer.size(); j++) {
        tempStacks[i][offset + j] = static_cast<OUTREAL>(tempBuffer[j]);
      }
    }
  }
  for (size_t i = 0; i < NXc; ++i) {
    double x = xmin + i * dh;
    for (size_t j = 0; j < NYc; ++j) {
      double y = ymin + j * dh;
      for (size_t k = 0; k < NZc; ++k) {
        double z = zmin + k * dh;
        size_t idx1d = k * NYc * NXc + j * NXc + i; //  x, y, z
        size_t idx3d = i * NYc * NZc + j * NZc + k; //  z, y, x
        Stacks[0][idx1d] = x;
        Stacks[1][idx1d] = y;
        Stacks[2][idx1d] = z;
        Stacks[3][idx1d] = tempStacks[0][idx3d];
        Stacks[4][idx1d] = tempStacks[1][idx3d];
        Stacks[5][idx1d] = tempStacks[2][idx3d];
      }
    }
  }
  std::cout << std::endl;
  //   if (showIb) {
  //     datasize = std::ceil(ndata / 4);
  //     std::vector<char> ibdata(datasize);
  //     infile.read((char *)&ibdata[0], datasize);
  //     ExtractIb(ibdata, &Stacks[7][offset], ndata);
  //   }
  //   infile.close();
  //   for (int i = zc_start; i <= zc_end; ++i)
  //     zids.erase(i);

  // if (!zids.empty()) {
  //   std::cout << "Inconsistent data size " << std::endl;
  //   std::cout << "missling z-slices: ";
  //   for (auto i : zids) {
  //     std::cout << i << " ";
  //   }
  //   std::cout << "." << std::endl;
  // }
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Error, insufficient input parameters, Usage: Combine "
                 "DatFlow/Flow0000.00000 0.plt"
              << std::endl;
    exit(-1);
  }
  // std::string optIblankStr("--iblank");
  // bool showIb = false;
  // for (int i = 1; i < argc; ++i) {
  //   if (optIblankStr.compare(argv[i]) == 0) {
  //     showIb = true;
  //     break;
  //   }
  // }
  // std::string optErrorStr("--error");
  // bool showError = false;
  // for (int i = 1; i < argc; ++i) {
  //   if (optErrorStr.compare(argv[i]) == 0) {
  //     showError = true;
  //     break;
  //   }
  // }
  std::string flow_file(argv[1]);
  std::string plt_file(argv[2]);
  std::vector<size_t> rawN(3);
  std::vector<std::vector<OUTREAL>> Stacks;
  std::vector<std::string> var = {"x", "y", "z", "u", "v", "w"};
  // if (showIb) {
  //   var.push_back("Ib");
  // }
  stackDataFromFiles(flow_file, Stacks, rawN[0], rawN[1], rawN[2]);
  OutputTec360_binary(plt_file, var, rawN, Stacks);
  std::cout << "Write " << plt_file << std::endl;
  // if (showError) {
  //   printErrors(mesh_file, var, Stacks);
  // }
  return 0;
}
