#include<iostream>
#include"F:/ENG/MONO ROBO/P/Net.h" /**modify to Net.h file path**/
#include<time.h>    //clock(),CLOCKS_PER_SEC for calculating execution time
frame readImage2(const char* filename);
Net* segmentation(frame F, Net* DN);

int main()
{
    clock_t start, end;    //for calculating Execution time
    frame F = readImage2("frame1.txt");
    Net* DN = new Net(F);
    start = clock();
    segmentation(F, DN);
    end=clock();
    std::cout<<"  "<<"Execution Time: "<<double(end-start)/double(CLOCKS_PER_SEC)<<" sec\n";
    DN->print_polys();
}
