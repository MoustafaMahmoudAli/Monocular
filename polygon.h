#ifndef _POLYGON_H
#define _POLYGON_H

#include<list>
#include"vertex.h"
class vertex;
#include"point.h"
#include"vector.h"
#include"frame.h"

class polygon{
    //private:
    public:
        //first list contains the vertices of the outer polygon in clockwise order
        //following lists -if any- contain child polygons in anti-clockwise order
        /*              ______________
                       |    _____     |
                       |   |     |    |
        outer polygon->|   |_____|    |
                       |       ^      |
                       |       |      |
                       |child polygon |
                       |______________|

        */
        std::list<std::list<vertex*>> myVertices;
        //number of vertices in the polygon
        int n;
        //average color of the polygon
        color c;
        //area covered by the polygon
        double area;

    //public:
        polygon();
        //construct a polygon from a list of vertices in clockwise order
        //first list contains the vertices' of the outer polygon in clockwise order
        //following lists -if any- contain child polygons in anti-clockwise order
        polygon(std::list<std::list<vertex*>>& L);
        //copy constructor
        polygon(const polygon& poly);
        //assignment operator
        polygon& operator=(const polygon& poly);
        bool operator == (const polygon& poly);
        //find next vertex after v in the clockwise order
        vertex* next_vertex(const vertex* v);
        //find previous vertex before v in the clockwise order
        vertex* prev_vertex(const vertex* v);
        //add vertex v to the polygon before vertex w
        void add_vertex(vertex* v,vertex* w);
        //remove vertex from the polygon
        void remove_vertex(const vertex* v);
        void add_poly_list(std::list<vertex*>& p);
        void remove_poly_list(std::list<vertex*>& p);
        std::list<std::list<vertex*>> list_of_list();

        //is the point in/on the polygon or outside
        //point in polygon
        bool PIP(point p);

        //calculate color average of the polygon
        color color_avg(frame F);

        //calculate the area of the polygon
        //Using Gauss's area formula (Shoelace formula)
        double compute_area();

    friend class Net;
};

#endif // _POLYGON_H

