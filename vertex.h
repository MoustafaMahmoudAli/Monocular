#ifndef _VERTEX_H_
#define _VERTEX_H_

#include<list>
#include"polygon.h"
class polygon;
#include"vector.h"
#include"point.h"

// a net vertex
class vertex{
    //private:
    public:
        //screen coordinates
        //origin at upper left corner
        //x increases to the right
        //y increases downward
        int x,y;
        //total deformation force;
        vector total_force;

        //list of pointers to polygons in which the vertex is a member in clockwise order
        std::list<polygon*> myPolygons;

    //public:
        vertex();
        //initialize vertex
        vertex(int x_coord, int y_coord);
        vertex(int x_coord, int y_coord, std::list<polygon*> polys);
        //copy vertex
        vertex(const vertex& v);
        //assignment operator
        vertex& operator=(const vertex& v);
        bool operator == (const vertex& v);
        //update x , y coordinates of the vertex
        void update_loc(int x_coord, int y_coord);
        //update the vertex's polygon list
        void update_polylist(std::list<polygon*>& polys);
        //insert a polygon p1 in myPolygon list before p2
        void add_polygon(polygon* p1, polygon* p2);
        //remove a polygon from myPolygon list
        void remove_polygon(polygon* p);
        // Return Point of The coordinates of The Vertex (x,y)
        point V_point();

    friend class polygon;
    friend class Net;
};

#endif // _VERTEX_H
