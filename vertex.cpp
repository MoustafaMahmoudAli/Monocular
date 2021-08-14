#include"vertex.h"

    vertex::vertex(){}

    //initialize vertex
    vertex::vertex(int x_coord, int y_coord) : x(x_coord), y(y_coord), total_force(){}
    vertex::vertex(int x_coord, int y_coord, std::list<polygon*> polys) : x(x_coord), y(y_coord), total_force(){
        myPolygons = polys;
    }

    //copy vertex
    vertex::vertex(const vertex& v){
        x = v.x;
        y = v.y;
        total_force = v.total_force;
        myPolygons = v.myPolygons;
    }
    //assignment operator
    vertex& vertex::operator=(const vertex& v){
        x = v.x;
        y = v.y;
        total_force = v.total_force;
        myPolygons = v.myPolygons;
        return *this;
    }

    bool vertex::operator == (const vertex& v){
        return ( (x==v.x) && (y==v.y) && (total_force==v.total_force) && (myPolygons == v.myPolygons) );
    }
    //update x , y coordinates of the vertex
    void vertex::update_loc(int x_coord, int y_coord){
        x = x_coord;
        y = y_coord;
    }
    //update the vertex's polygon list
    void vertex::update_polylist(std::list<polygon*>& polys){
        myPolygons = polys;
    }
    //insert a polygon p1 in myPolygon list before p2
    void vertex::add_polygon(polygon* p1, polygon* p2){
        typedef std::list<polygon*>::iterator poly_iter;
        for(poly_iter i = myPolygons.begin(); i != myPolygons.end(); ++i){
            if(*i == p2){
                myPolygons.insert(i, p1);
                return;
            }
        }
    }
    //remove a polygon from myPolygon list
    void vertex::remove_polygon(polygon* p){
        myPolygons.remove(p);
    }

    // Return Point of The coordinates of The Vertex (x,y)
    point vertex::V_point()
    {
        return point( x, y);
    }
