#ifndef _DEF_NET_H_
#define _DEF_NET_H_

#include<list>
#include"vertex.h"
#include"polygon.h"
#include"point.h"
#include"vector.h"
#include"frame.h"
#include<iostream>
#include<fstream>

    //the two disparities associated with vertex v2
    struct ver_disp{
        //disparity for vertex v2 with the next vertex v3
        //   ______
        //  |      |
        //v2--------------->v3
        int H_next;
        //disparity for vertex v2 with the previous vertex v1
        //           ______
        //          |      |
        //v1--------------->v2
        int H_prev;
    };

    //the four disparities associated with an edge v1--->v2
    struct edge_disp{
        //right disparity for this polygon
        //           ______
        //          |      |
        //v1--------------->v2
        int H_R1;
        //left disparity for this polygon
        //   ______
        //  |      |
        //v1--------------->v2
        int H_L1;
        //right disparity for the neighboring polygon
        //v1--------------->v2
        //  |______|
        int H_R2;
        //left disparity for the neighboring polygon
        //v1--------------->v2
        //          |______|
        int H_L2;
    };

    // Sensitivity Region diplacement Values
    struct diplacement
    {
        double deltaX;
        double deltaY;
        diplacement():deltaX(0.0),deltaY(0.0){}
    };

    struct scaning_d
    {
        int X_max;
        int X_min;
        int Y_max;
        int Y_min;

    };

    struct Edge
    {
        vertex *S;   // start
        vertex *E;   //End
    };


//Deformable Net
class Net{
    //private:
    public:
        //list of all vertices in the net
        std::list<vertex*> verts;
        //list of all polygons in the net
        std::list<polygon*> polys;
        //list of vertices that has not finished deformation yet
        std::list<vertex*> deformable;
        //image dimensions
        int image_width, image_height;
        /*threshold parameters*/
        int maintenance_period;
        //sensitivity region width
        int S_width;
        //color threshold
        int zeta_thresh;
        //merge size threshold
        int size_thresh;
        //overall small force threshold (for vertices insertion)
        //this parameter is tuned experimentally
        float force_thresh;
        //vertex insertion threshold(percentage of sensitivity region area)
        float segma_thresh;
        //filling factor threshold
        float fill_thresh;

    //public:
        //Net initialization to image dimensions (width x hight)
        //assuming 10x10 pixels squares (last column and last row may have smaller polygons)
        Net(frame F);
        //copy constructor
        Net(const Net& net_obj);
        //assignment operator
        Net& operator=(const Net& net_obj);

        //ceil function
        int ceil(float N);

        void add_polygon(polygon *p);
        void remove_polygon(polygon *p);
        void remove_vertex(vertex* v);

        //=====================================================================================
        /**Deformation Cycle**/

        //get neighbor of a polygon with edge v1-->v2
        //v1 precede v2 in the clockwise order
        polygon* get_nbour(polygon* poly, vertex* v1, vertex* v2);

        //calculate the color distance between a polygon average color and a pixel color
        //cylindrical distance
        double dist_cyl(color i, color j);
        //chromaticity distance (dist_cyl help function)
        double dist_chrom(color i,color j);

        //calculate the disparity for vertex v2 with both next vertex v3 and previous vertex v1 in the polygon
        ver_disp vertex_disparity(polygon* poly, vertex* v1, vertex* v2, vertex* v3,frame F);
        //(vertex_disparity help funcitons)
        diplacement measure_diplacement(vertex *v1, vertex *v2);
        scaning_d measure_scaning_D (vertex *v1 ,vertex *v2 ,vertex *v3 ,vertex *v4);
        void coordinates_check(int* x, int* y, int max_X, int max_Y);

        //go through all vertices in all polygons
        //calculate the deformation forces and update the location of the vertices
        //completing a full deformation cycle
        //return the total energy in a single deformation cycle
        //Energy is the number of vertices that changed place in one deformation cycle
        float deformation(frame F);

        //==============================================================================================
        /**Maintenance Cycle**/

        //merge two neighboring polygons
        //return a pointer to the merged polygon
        //updates all affected vertices and polygons in the Net
        polygon* merge_utility(polygon* p1,polygon* p2, frame F);
        //(merge_utility help function)
        bool inList(std::list<vertex*> L, vertex* v);

        //merge neighboring polygons having average color difference below color threshold
        void color_merge(frame F);
        //merge a polygon with an area <= size threshold to a neighboring polygon having the nearest average color
        void size_merge(frame F);

        //split contracted polygons
        void split_polygons(frame F);
        //(split_polygons help functions)
        bool connected(polygon *P,Edge *e1, Edge*e2);
        bool child_Polygon_Exist(std::list<std::list<vertex*>> vertices, polygon* P);
        bool edges_contracted(Edge *e1, Edge *e2,frame&F);

        //calculate the filling factor for each polygon in the net
        //A polygon having a filling factor below 80% is reinitialized
        void insert_polygons(frame F);

        //delete unnecessary vertices
        void delete_vertices();

        //break edges that have an overall force below the force threshold with a disparity measure above segma threshold
        void insert_vertices(frame F);
        //(insert_vertices help function)
        //calculate the four disparities associated with an edge v1--->v2 in poly
        edge_disp edge_disparity(polygon* poly, polygon* nbour, vertex* v1, vertex* v2, frame F);

        //A full maintenance cycle
        void maintenance(frame F);

        /**Debugging**/
        void print_polys();
        void print_verts();

};

#endif // _DEF_NET_H
