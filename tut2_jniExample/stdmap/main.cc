#include <map>
#include <cstdio>

int main(int argc,char** argv){
    std::map<int*,int> testMap;
    int tmp1=123;
    int tmp2=1123;
    int tmp3=2123;
    testMap[&tmp1]=456;
    testMap[&tmp2]=56;
    testMap[&tmp3]=6;
    return 0;
}
