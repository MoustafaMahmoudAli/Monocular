#include"point.h"

    point::point():x(0), y(0) {}

    point::point(int x_coord, int y_coord): x(x_coord), y(y_coord) {}

    point::point(const point& p){
        x = p.x;
        y = p.y;
    }

    point& point::operator = (const point& p){
        x = p.x;
        y = p.y;
        return *this;
    }

    bool point::operator ==(const point& p){
        return (x==p.x && y==p.y);
    }

    bool point::operator !=(const point& p){
        return (x!=p.x || y!=p.y);
    }

    void point::update_coords(int x_coord, int y_coord){
        x = x_coord;
        y = y_coord;
    }

    point point::move(vector vec){
        return point(x+vec.x, y+vec.y);
    };
