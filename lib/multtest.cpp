/*
#include "mslib.h"
#include "testlib.h"
*/
#include <chrono>
#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char**  argv) {
    cout<<"Hello World!";
    auto t1=chrono::high_resolution_clock::now();
    //vector<struct SqMat> testMat=genRandUnitValElmSqMatList(10,50);
    auto t2=chrono::high_resolution_clock::now();
    cout<< " "
              << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n"<<endl;
    return 1;
}
