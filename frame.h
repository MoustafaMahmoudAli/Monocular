#ifndef _FRAME_H
#define _FRAME_H

/*
by Shady R. Fam.
*/

//HSI Color
class color
{
    //private:
    public:
        float R, G, B;
        float H, S, I;
    //public:
        color();
        color(float r, float g, float b);
        color(const color& c);
        color& operator = (const color& c);
        bool operator == (const color& c);
    friend class Net;
};

//Image frame, only contains dimensions and HSI pixel data
struct frame
{
        color** data;
        int height,width;
};

#endif // _FRAME_H
