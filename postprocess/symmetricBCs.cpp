#include<iostream>
#include<vector>
#include<string>
#include<cstdio>
using namespace std;

vector<int> dx{0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0};
vector<int> dy{0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1};
vector<int> dz{0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1};

int findid(int x, int y, int z)
{
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(x==dx[i] && y==dy[i] && z==dz[i])
        {
            return i;
        }
    }
    return -1;
}

int main()
{
    string format("!%3d%3d%3d%3d -> %3d%3d%3d%3d\n");
    // wall boundary
    cout << "x min wall" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dx[i]==-1)
        {
            int x = -dx[i];
            int y = -dy[i];
            int z = -dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    cout << "x max wall" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dx[i]==1)
        {
            int x = -dx[i];
            int y = -dy[i];
            int z = -dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    cout << "y min wall" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dy[i]==-1)
        {
            int x = -dx[i];
            int y = -dy[i];
            int z = -dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    cout << "y max wall" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dy[i]==1)
        {
            int x = -dx[i];
            int y = -dy[i];
            int z = -dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    cout << "z min wall" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dz[i]==-1)
        {
            int x = -dx[i];
            int y = -dy[i];
            int z = -dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    cout << "z max wall" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dz[i]==1)
        {
            int x = -dx[i];
            int y = -dy[i];
            int z = -dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    // symmetric boundary
    cout << "x min symmetric" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dx[i]==-1)
        {
            int x = -dx[i];
            int y = dy[i];
            int z = dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    cout << "x max symmetric" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dx[i]==1)
        {
            int x = -dx[i];
            int y = dy[i];
            int z = dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    cout << "y min symmetric" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dy[i]==-1)
        {
            int x = dx[i];
            int y = -dy[i];
            int z = dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    cout << "y max symmetric" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dy[i]==1)
        {
            int x = dx[i];
            int y = -dy[i];
            int z = dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    cout << "z min symmetric" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dz[i]==-1)
        {
            int x = dx[i];
            int y = dy[i];
            int z = -dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
    cout << "z max symmetric" << endl;
    for(size_t i=0; i<dx.size(); ++i)
    {
        if(dz[i]==1)
        {
            int x = dx[i];
            int y = dy[i];
            int z = -dz[i];
            printf(format.c_str(), i, dx[i], dy[i], dz[i], findid(x, y, z), x, y, z);
        }
    }
}
