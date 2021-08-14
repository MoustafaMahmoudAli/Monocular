#include"Net.h"
#include<time.h>

    /**automatic color image segmentation**/
    //returns pointer to the segmented Net
    //Prepared BY : Khaled M. Kamel
    void segmentation(frame F, Net* DN){
        typedef std::list<vertex*>::iterator ver_iter;
        //Merge the neighboring polygons that have similar average colors
        DN->color_merge(F);
        //DN->size_merge(F);
        //remove all vertices on the image border that are part of only a single polygon except the corners
        vertex* v;
        int x,y;
        std::list<vertex*> temp = DN->verts;
        for(ver_iter i = temp.begin(); i != temp.end(); ++i){
            v = *i;
            x = v->x;
            y = v->y;
            if( v->myPolygons.size() == 1 &&
               !( (x==0&&y==0)||(x==0&&y==F.height)||(x==F.width&&y==0)||(x==F.width&&y==F.height) )){
                v->myPolygons.front()->remove_vertex(v);
                DN->remove_vertex(v);
            }
        }
       //maximum vertex displacement in one deformation cycle
        float disp;
        //number of deformation cycles
        int iterations = 1;
            while(true){
                disp = DN->deformation(F);
                if(disp > 0){
                    if(iterations % DN->maintenance_period == 0){
                        DN->maintenance(F);
                    }
                    ++iterations;
                }
                else
                    break;
            }
            DN->delete_vertices();
    }
