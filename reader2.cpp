#include <fstream>
#include"F:/ENG/MONO ROBO/P/frame.h"
#include<iostream>
#include<cmath>
#include<time.h>
    void computeHSI(color& c){
        float r,g,b;
        r=c.R; r/=255;
        g=c.G; g/=255;
        b=c.B; b/=255;
        c.I = (r+g+b)/3;
        c.I *= 100;
        float pmin = std::min(r,std::min(g,b));
        if(r == g && g == b)
            c.S = 0;
        else
            c.S = (c.I == 0) ? 0 : 1 - (3/(r+g+b)*pmin);
        c.S *= 100;
        if(c.S == 0)
            c.H = 0;
        else{
            c.H = .5*(2*r-g-b) / sqrt((r-g)*(r-g) + (r-b)*(g-b));
            if(c.H > 1)
                c.H = 1;
            else if (c.H < -1)
                c.H = -1;
            c.H = acos(c.H) * 57.296;
        }
    }

frame readImage2(const char* filename){

    frame F;
    std::ifstream file(filename);
    file >> F.width >> F.height;
    F.data = new color*[F.width];
    for (int i = 0; i < F.width; i++){
        F.data[i] = new color[F.height];
    }
    int R, G, B;
    for(int j = 0; j < F.height; ++j){
        for(int i = 0; i < F.width; ++i){
            file>>R>>G>>B;
            F.data[i][j].R = R;
            F.data[i][j].G = G;
            F.data[i][j].B = B;
            computeHSI(F.data[i][j]);
        }
    }
    return F;
}
