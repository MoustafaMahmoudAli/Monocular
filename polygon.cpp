#include"polygon.h"
#include<iostream>
#include<fstream>
#include"point.h"
#include"vertex.h"
#include<cmath>/**is this allowed?**/
#include<time.h>

void computeHSI(color& c);

    polygon::polygon(): n(0), c(),area(0) {}

    //construct a polygon from a list of vertices in clockwise order
    //first list contains the vertices' of the outer polygon in clockwise order
    //following lists -if any- contain child polygons in anti-clockwise order
    polygon::polygon(std::list<std::list<vertex*>>& L): c(),area(0){
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        myVertices = L;
        n = 0;
        for(list_iter i = L.begin(); i != L.end(); ++i){
            n += (*i).size();
        }
    }

    //copy constructor
    polygon::polygon(const polygon& poly){
        myVertices = poly.myVertices;
        n = poly.n;
        c = poly.c;
        area = poly.area;
    }

    //assignment operator
    polygon& polygon::operator=(const polygon& poly){
        myVertices = poly.myVertices;
        n = poly.n;
        c = poly.c;
        area = poly.area;
        return *this;
    }

    bool polygon::operator == (const polygon& poly){
        return ( (myVertices==poly.myVertices) && (n==poly.n) && (c==poly.c) && (area==poly.area) );
    }

    //find next vertex after v
    vertex* polygon::next_vertex(const vertex* v){
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;
        //search for vertex v
        for(list_iter i = myVertices.begin(); i != myVertices.end(); ++i){
            for(ver_iter j = (*i).begin(); j != (*i).end(); ++j){

                if(*j == v){//return the next vertex
                    return (++j == (*i).end()) ? *((*i).begin()) : *j;
                }
            }
        }
    }

    //find previous vertex before v
    vertex* polygon::prev_vertex(const vertex* v){
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;
        //search for vertex v
        for(list_iter i = myVertices.begin(); i != myVertices.end(); ++i){
            for(ver_iter j = (*i).begin(); j != (*i).end(); ++j){

                if(*j == v){//return the previous vertex
                    return (--j == (*i).end()) ? *(--(*i).end()) : *j;
                }
            }
        }
    }

    //add vertex v to the polygon before vertex w
    void polygon::add_vertex(vertex* v,vertex* w){
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;
        //search for vertex w
        for(list_iter i = myVertices.begin(); i != myVertices.end(); ++i){
            for(ver_iter j = (*i).begin(); j != (*i).end(); ++j){

                if(*j == w){//add vertex v before vertex w
                    (*i).insert(j,v);
                    ++n;
                    return;
                }
            }
        }
        std::cout<<"Add before a non-existent vertex";
    }

    //remove vertex from the polygon
    void polygon::remove_vertex(const vertex* v){
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;
        //search for vertex v
        for(list_iter i = myVertices.begin(); i != myVertices.end(); ++i){
            for(ver_iter j = (*i).begin(); j != (*i).end(); ++j){

                if(*j == v){//remove vertex v
                    (*i).erase(j);
                    --n;
                    return;
                }
            }
        }
//        std::cout<<"Removal of a non-existent vertex";
    }

    void polygon::add_poly_list( std::list<vertex*>& p)
    {
      myVertices.push_back( p);
    }

    void polygon::remove_poly_list(  std::list<vertex*>& p )
    {
       myVertices.remove( p);
    }

    std::list<std::list<vertex*>> polygon:: list_of_list()
    {
        return myVertices;
    }

    //Prepared BY : Ahmed Maged
    bool polygon::PIP(point p)
    {
        // Test the ray against all sides
        int intersections = 0;
        std::list<vertex*>::iterator next;
        for (auto& j : myVertices)
        {
            next = j.begin();
            next++;

            int  maxY = INT_MIN, minY = INT_MAX, currentY = INT_MAX;
            for (auto& k : j)
            {
                if (k->y > maxY)
                    maxY = k->y;
                if (k->y < minY)
                    minY = k->y;
            }
            for (auto i = j.begin(); i != j.end(); i++, next++)
            {
                if (next == j.end())
                {
                    next = j.begin();
                    if ((*next)->y == p.y)
                        continue;
                }
                // Test if current side intersects with ray.
                // tes y cor.
                if ((p.y >= (*i)->y && p.y <= (*next)->y) || (p.y <= (*i)->y && p.y >= (*next)->y))
                {
                    //on the vertex
                    if (((*i)->x == p.x && (*i)->y == p.y) || ((*next)->x == p.x && (*next)->y == p.y))
                    {
                        return 1;
                    }
                    //point on right
                    if ((*i)->x < p.x && (*next)->x < p.x)
                    {
                        continue;
                    }

                    //point on left
                    if ((*i)->x > p.x && (*next)->x > p.x)
                    {
                        if (p.y == maxY || p.y == minY)
                            continue;
                        if ((*i)->y == currentY && p.y == currentY)
                            continue;

                        currentY = (*next)->y;

                        intersections++;
                        continue;
                    }

                    else //between
                    {
                        if ((*next)->y == (*i)->y)
                            return 1;
                        if ((*next)->x == (*i)->x)
                            return 1;
                        if (p.y == maxY || p.y == minY)
                            continue;
                        if ((*i)->y == currentY && p.y == currentY)
                            continue;
                        currentY = (*next)->y;
                        float m = ((*next)->y - (*i)->y) / (float)((*next)->x - (*i)->x);
                        float c = (*next)->y - (m * (*next)->x);
                        float xpos = (p.y - c) / m;
                        if (xpos > p.x)
                            intersections++;
                        else if (xpos == p.x) //on the edge
                            return 1;

                    }
                }
                else
                {
                    currentY = INT_MAX;
                }
            }

        }
        if ((intersections & 1) == 1)
            return 1;
        else
            return 0;
    }

    //BY: Shay R. Fam
    color polygon::color_avg(frame F)
    {

        point temp;
        std::list<vertex*> L=*(myVertices.begin());
        int avgR=0; int avgG=0; int avgB=0; int pixelCount=0;
        int minx=F.width; int maxx=0; int miny=F.height; int maxy=0;
        std::list<vertex*>::iterator it;
        for (it=L.begin();it!=L.end();++it)
        {
            if ((*it)->x < minx) minx=(*it)->x;
            if ((*it)->x > maxx) maxx=(*it)->x;
            if ((*it)->y < miny) miny=(*it)->y;
            if ((*it)->y > maxy) maxy=(*it)->y;
        }

        for (int i=minx;i<=maxx;i++)
        {
            temp.x=i;
            for (int j=miny;j<=maxy;j++)
            {
                temp.y=j;
                if (polygon::PIP(temp))
                {
                    avgR+=F.data[i][j].R;
                    avgG+=F.data[i][j].G;
                    avgB+=F.data[i][j].B;
                    pixelCount++;
                }
            }
        }
        avgR/=pixelCount; avgG/=pixelCount; avgB/=pixelCount;
        color c;
        c.R=avgR; c.G=avgG; c.B=avgB;
        computeHSI(c);

        return c;
    }

    //return a pixel having the average color of the polygon
    //Prepared BY: Khaled M. Kamel
//    pixel polygon::RGBcolor(){
////        /**debug**/
////        std::cout<<c.H<<" "<<c.S<<" "<<c.I<<'\n';
////        /**debug**/
////        float H = c.H/60;
////        float I = c.I/100;
////        float S = c.S/100;
////        float Z = 1 - std::abs( (int(H)%2) - 1 );
////        float C = (3*I*S)/(1+Z);
////        float X = C*Z;
////        float R, G, B;
////        if(H >=0 && H <= 1){
////            R = C;
////            G = X;
////            B = 0;
////        }
////        else if(H >=1 && H <= 2){
////            R = X;
////            G = C;
////            B = 0;
////        }
////        else if(H >=2 && H <= 3){
////            R = 0;
////            G = C;
////            B = X;
////        }
////        else if(H >=3 && H <= 4){
////            R = 0;
////            G = X;
////            B = C;
////        }
////        else if(H >=4 && H <= 5){
////            R = X;
////            G = 0;
////            B = C;
////        }
////        else if(H >=5 && H <= 6){
////            R = C;
////            G = 0;
////            B = X;
////        }
////        float m = I*(1-S);
////        pixel p( 255*(R+m), 255*(G+m), 255*(B+m) );
////
//////        /**debug**/
//////        if(p.R < 0 || p.R > 255 || p.G < 0 || p.G > 255 || p.B < 0 || p.B > 255){
//////            int s = 0;
//////        }
//////        /**debug**/
//        float H = c.H, S = c.S/100, I = c.I/100;
//        float R, G, B;
//        if(H >= 0 && H <= 120){
//            B = I*(1-S);
//            R = I*(1 + (S*cos(H/57.296))/cos((60-H)/57.296));
//            G = 3*I - R - B;
//        }
//        else if(H > 120 && H <= 240){
//            H -= 120;
//            R = I*(1-S);
//            G = I*(1 + (S*cos(H/57.296))/cos((60-H)/57.296));
//            B = 3*I - R - G;
//        }
//        else if(H > 240 && H <= 360){
//            H -= 240;
//            G = I*(1-S);
//            B = I*(1 + (S*cos(H/57.296))/cos((60-H)/57.296));
//            R = 3*I - G - B;
//        }
//        pixel p(255*R, 255*G, 255*B);
////        /**debug**/
////        std::cout<<p.R<<" "<<p.G<<" "<<p.B<<'\n';
////        int s = 0;
////        /**debug**/
//
//        return p;
//    }

    //calculate the area of the polygon
    //Using Gauss's area formula (Shoelace formula)
    //Prepared BY: Khaled M. Kamel
    double polygon::compute_area(){
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;
        //Area of polygon i is the area of the outer polygon subtracting the areas of all child polygons
        //number of partial areas
        int n = myVertices.size();
        //initialize partial areas
        double part_area[n] = {0};
        vertex *v1, *v2;
        int current_area = 0;
        //calculate partial areas
        for(list_iter i = myVertices.begin(); i != myVertices.end(); ++i){
            for(ver_iter j = (*i).begin(); j != (*i).end(); ++j){
                //current vertex
                v1 = *j;
                //next vertex
                v2 = (++j == (*i).end()) ? *(*i).begin() : *j;
                --j;
                part_area[current_area] += v1->x*v2->y - v1->y*v2->x;
            }
            part_area[current_area] = 0.5 * std::abs(part_area[current_area]);
            ++current_area;
        }
        //calculate total polygon area
        for(int i = 1; i < n; ++i){
            part_area[0] -= part_area[i];
        }
        return part_area[0];
    }
