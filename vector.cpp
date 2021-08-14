#include"vector.h"
#include<math.h>
    vector::vector() : x(0), y(0) {}

    vector::vector(float x_coord,float y_coord) : x(x_coord), y(y_coord) {}

    vector::vector(const vector& vec){
        x = vec.x;
        y = vec.y;
    }

    vector& vector::operator=(const vector& vec){
        x = vec.x;
        y = vec.y;
        return *this;
    }

    bool vector::operator == (const vector& vec){
        return ( (x==vec.x) && (y==vec.y) );
    }

    void vector::update_coords(float x_coord, float y_coord){
        x = x_coord;
        y = y_coord;
    }

    vector vector::operator + (vector vec) const{
        return vector(x+vec.x, y+vec.y);
    }

    vector vector::operator - (vector& vec) const{
        return vector(x-vec.x, y-vec.y);
    }

    void vector::operator += (vector& vec){
        x += vec.x;
        y += vec.y;
    }

    void vector::operator -= (vector& vec){
        x -= vec.x;
        y -= vec.y;
    }

    vector vector::operator * (float scale) const{
        return vector(x*scale, y*scale);
    }

    float vector::mag() const{
        return sqrt(x*x+y*y);
    }

    vector vector::unit_vec() const{
        return vector(x/mag(), y/mag());
    }

    vector vector::perpendicular() const{
        return vector(y,-x);
    }

    //dot product
    float vector::dot(vector& vec) const{
        return (x*vec.x + y*vec.y);
    }
