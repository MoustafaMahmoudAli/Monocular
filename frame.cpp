#include"frame.h"
#include<iostream>
#include<cmath>

#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)

    color::color() : R(0), G(0), B(0), H(0), S(0), I(0){}

    color::color(float r, float g, float b):R(r), G(g), B(b) {}

    color::color(const color& c){
        R = c.R;
        G = c.G;
        B = c.B;
        H = c.H;
        S = c.S;
        I = c.I;
    }

    color& color::operator = (const color& c){
        R = c.R;
        G = c.G;
        B = c.B;
        H = c.H;
        S = c.S;
        I = c.I;
        return *this;
    }

    bool color::operator == (const color& c){
        return (H==c.H && S==c.S && I==c.I);
    }


