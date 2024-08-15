#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

int main(int argc, char* argv[])
{
    string filename(argv[1]);
    string ofilename(argv[2]);
    ifstream ifile(filename.c_str());
    ofstream ofile(ofilename.c_str());
    char buffer[1000];
    ifile.getline(buffer,sizeof(buffer));
    string match("(iFish,");
    while(!ifile.eof())
    {
        string line(buffer);
        while(1)
        {
            size_t pos = line.find(match);
            if(pos!=line.npos)
            {
                size_t i, count = 1;
                for(i=pos+1; i<line.size();++i)
                {
                    if(line[i]=='(') ++count;
                    else if(line[i]==')') --count;
                    if(count == 0) break;
                }
                for(size_t j=pos+7; j<i; ++j)
                {
                    line[j-6] = line[j];
                }
                line[i-6] = ',';
                line[i-5] = 'i';
                line[i-4] = 'F';
                line[i-3] = 'i';
                line[i-2] = 's';
                line[i-1] = 'h';
            }
            else break;
        }
        ofile << line << endl;
        ifile.getline(buffer,sizeof(buffer));
    }
    ofile << buffer << endl;
}
