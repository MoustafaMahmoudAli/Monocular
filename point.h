#ifndef _POINT_H
#define _POINT_H
#include"vector.h"

class point{
    //private:
    public:
        int x, y;
    //public:
        point();
        point(int x_coord, int y_coord);
        point(const point& p);
        point& operator = (const point& p);
        bool operator ==(const point& p);
        bool operator !=(const point& p);
        void update_coords(int x_coord, int y_coord);
        point move(vector vec);

    friend class polygon;
    friend class Net;
};

#endif // _POINT_H
