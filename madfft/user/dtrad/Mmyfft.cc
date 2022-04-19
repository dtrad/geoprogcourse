#include "rsf.hh"
#include <vector>
#include <iostream>
#include <string>
#include "ffthelper.hh"

using namespace std;

int main(int argc, char** argv) {
    sf_init(argc, argv); // Initialize RSF

    iRSF par(0);
    vector<float> bandpass(4);    
    
    par.get("lcut",bandpass[0],1.);
    par.get("lpass",bandpass[1],2.);
    par.get("hpass",bandpass[2],10.);
    par.get("hcut",bandpass[3],14.);
    
    // Get input
    iRSF input;
    int n1, n2, n3;
    float d1, d2, d3;

    input.get("n1", n1);
    input.get("n2", n2);
    input.get("n3", n3);
    input.get("d1", d1);
    input.get("d2", d2);
    input.get("d3", d3);
    

    //string title= sf_getstring("title");
    // Set output
    oRSF output;
    std::valarray<float> in(n1 * n2);
    std::valarray<float> out(n1 * n2);
    output.put("n3", n3);
    int size = in.size();
    
    FFThelper* fft=new FFThelper(n1,n2,d1,bandpass);
    int ishot = 0;
    for (int i3 = 0; i3 < n3; i3++) {
        // find location for the next shot        
        input.seek(i3 * size * sizeof (float), SEEK_SET);
        input >> in;
        // fft each trace
        fft->apply(in,n2,1);
        fft->filter();
        fft->apply(out,n2,0);
        //out=in;
        cerr << "shot=" << ishot << " ";
        output << out;
        ishot++;
    }
    
    delete fft;
    exit(0);
}

