#ifndef _VECTOR_H
#define _VECTOR_H

class vector{
    //private:
    public:
        float x, y;
    //public:
        vector();
        vector(float x_coord,float y_coord);
        vector(const vector& vec);
        vector& operator=(const vector& vec);
        bool operator == (const vector& vec);
        void update_coords(float x_coord, float y_coord);
        vector operator + (vector vec) const;
        vector operator - (vector& vec) const;
        void operator += (vector& vec);
        void operator -= (vector& vec);
        vector operator * (float scale) const;
        //return the magnitude of the vector
        float mag() const;
        //return a unit vector in direction of this vector
        vector unit_vec() const;
        //return a perpendicular to this vector (+90 vector)
        vector perpendicular() const;
        //dot product
        float dot(vector& vec) const;
    friend class Net;
};

#endif // _VECTOR_H
