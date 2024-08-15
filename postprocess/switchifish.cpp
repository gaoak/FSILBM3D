#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

string trimedstr(char buffer[])
{
    string empty(" \r\n\t");
    string line(buffer);
    size_t pos = line.find_last_not_of(empty);
    if(pos!=line.npos)
    {
        return line.substr(0,pos+1);
    }
    else
    {
        return string("");
    }
}

int main(int argc, char* argv[])
{
    string filename(argv[1]);
    string ofilename(argv[2]);
    ifstream ifile(filename.c_str());
    ofstream ofile(ofilename.c_str());
    char buffer[1000];
    ifile.getline(buffer,sizeof(buffer));
    string match1("(iFish,");
    string match2("(nFish,");
    while(!ifile.eof())
    {
        string line = trimedstr(buffer);
        while(1)
        {
            size_t pos = line.find(match1);
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
        while(1)
        {
            size_t pos = line.find(match2);
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
                line[i-5] = 'n';
                line[i-4] = 'F';
                line[i-3] = 'i';
                line[i-2] = 's';
                line[i-1] = 'h';
            }
            else break;
        }
        string match3("(1:nFish,");
        while(1)
        {
            size_t pos = line.find(match3);
            if(pos!=line.npos)
            {
                size_t i, count = 1;
                for(i=pos+1; i<line.size();++i)
                {
                    if(line[i]=='(') ++count;
                    else if(line[i]==')') --count;
                    if(count == 0) break;
                }
                for(size_t j=pos+9; j<i; ++j)
                {
                    line[j-8] = line[j];
                }
                line[i-8] = ',';
                line[i-7] = '1';
                line[i-6] = ':';
                line[i-5] = 'n';
                line[i-4] = 'F';
                line[i-3] = 'i';
                line[i-2] = 's';
                line[i-1] = 'h';
            }
            else break;
        }
        if(line.size())
        {
            ofile << line << endl;
        }
        else
        {
            ofile << endl;
        }
        ifile.getline(buffer,sizeof(buffer));
    }
    string line = trimedstr(buffer);
    if(line.size())
    {
        ofile << line << endl;
    }
}
