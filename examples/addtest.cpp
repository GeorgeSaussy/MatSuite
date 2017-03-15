#include "mslib.h"
#include "testlib.h"
#include <chrono>
#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char**  argv) {
    cout<<"Hello World!"<<endl;
    vector<struct SqMat> testMat=genRandUnitValElmSqMatList(50,10);
    int k=0;
    int k1=0;
    struct SqMat ongoing;
    initZeroSqMat(&ongoing,10);
    auto t1=chrono::high_resolution_clock::now();
    for(k=0;k<50;k++) {
        for(k1=k;k<50;k++) {
            addSafeSqMat(testMat[k],testMat[k1],&ongoing);
        }
    }
    auto t2=chrono::high_resolution_clock::now();
    cout<< " "
              << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n"<<endl;
    return 1;
}
