#include"Net.h"
#include<cmath> /**is this allowed?**/
#include<algorithm> /**TOTALY NOT ALLOWED**/
#include<time.h>
void computeHSI(color& c);

    //Net initialization to image dimensions (width x hight)
    //assuming 10x10 pixels squares (last column and last row may have smaller polygons)
    Net::Net(frame F):maintenance_period(20), S_width(10), zeta_thresh(20), size_thresh(150),
                       force_thresh(1.5), segma_thresh(0.2), fill_thresh(.8){

        image_width = F.width;
        image_height = F.height;
        /*generating vertices*/
        int x = 0, y = 0;
        //number of squares in x and y directions
        int x_squares = ceil(float(image_width)/10), y_squares = ceil(float(image_height)/10);
        //temp array used to fill polygons' vertex list
        vertex* temp_verts[y_squares+1][x_squares+1];
        //fill the vertices list
        for(int j = 0; j < y_squares; ++j){
            x = 0;
            for(int i = 0; i < x_squares; ++i){
                vertex* v = new vertex(x,y);
                verts.push_back(v);
                //save pointer to the vertex
                temp_verts[j][i] = v;
                x += 10;
            }
            //add last vertex in row
            vertex* v = new vertex(image_width,y);
            verts.push_back(v);
            temp_verts[j][x_squares] = v;
            y += 10;
        }
        //add last row
        x = 0;
        for(int i = 0; i < x_squares; ++i){
            vertex* v = new vertex(x,image_height);
            verts.push_back(v);
            temp_verts[y_squares][i] = v;
            x += 10;
        }
        vertex* v = new vertex(image_width,image_height);
        verts.push_back(v);
        temp_verts[y_squares][x_squares] = v;
        //temp array used to fill vertices' polygon list in clockwise order
        polygon* temp_polys[y_squares][x_squares];
        /*generating polygons*/
        int avgR, avgG, avgB;
        int h, w;
        for(int j = 0; j < y_squares; ++j){
            for(int i = 0; i < x_squares; ++i){
                //create new polygon
                std::list<std::list<vertex*>> L;
                std::list<vertex*> temp;
                L.push_front(temp);
                (*L.begin()).push_back(temp_verts[j][i]);
                (*L.begin()).push_back(temp_verts[j][i+1]);
                (*L.begin()).push_back(temp_verts[j+1][i+1]);
                (*L.begin()).push_back(temp_verts[j+1][i]);
                polygon* p = new polygon(L);
                polys.push_back(p);
                temp_polys[j][i] = p;
                if(j == y_squares-1 && i == x_squares-1){
                    h = (F.height%10);
                    w = (F.width%10);
                }
                else if(i == x_squares-1){
                    h = 10;
                    w = (F.width%10);
                }
                else if(j == y_squares-1){
                    h = (F.height%10);
                    w = 10;
                }
                else{
                    h = 10;
                    w = 10;
                }
                p->area = h * w;
                h += 10*j;
                w += 10*i;
                avgR = 0; avgG = 0; avgB = 0;
                for(int y = 10*j; y < h; ++y){
                    for(int x = 10*i; x < w; ++x){
                        avgR += F.data[x][y].R;
                        avgG += F.data[x][y].G;
                        avgB += F.data[x][y].B;
                    }
                }
                color avg(avgR/p->area, avgG/p->area, avgB/p->area);
                computeHSI(avg);
                p->c = avg;
            }
        }
        //fill vertices' polygon list in clockwise order
        for(int j = 0; j <= y_squares; ++j){
            for(int i = 0; i <= x_squares; ++i){
                //assuming vertex at the center
                //upper left polygon
                if(j != 0 && i != 0)
                    temp_verts[j][i]->myPolygons.push_back(temp_polys[j-1][i-1]);
                //upper right polygon
                if(j != 0 && i != x_squares)
                    temp_verts[j][i]->myPolygons.push_back(temp_polys[j-1][i]);
                //lower right polygon
                if(j != y_squares && i != x_squares)
                    temp_verts[j][i]->myPolygons.push_back(temp_polys[j][i]);
                //lower left polygon
                if(j != y_squares && i != 0)
                    temp_verts[j][i]->myPolygons.push_back(temp_polys[j][i-1]);
            }
        }
        deformable = verts;
    }

    //copy constructor
    Net::Net(const Net& net_obj){
        verts = net_obj.verts;
        polys = net_obj.polys;
        deformable = net_obj.deformable;
        image_width = net_obj.image_width;
        image_height = net_obj.image_height;
        maintenance_period = net_obj.maintenance_period;
        S_width = net_obj.S_width;
        zeta_thresh = net_obj.zeta_thresh;
        size_thresh = net_obj.size_thresh;
        force_thresh = net_obj.force_thresh;
        segma_thresh = net_obj.segma_thresh;
        fill_thresh = net_obj.fill_thresh;
    }

    //assignment operator
    Net& Net::operator=(const Net& net_obj){
        verts = net_obj.verts;
        polys = net_obj.polys;
        deformable = net_obj.deformable;
        image_width = net_obj.image_width;
        image_height = net_obj.image_height;
        maintenance_period = net_obj.maintenance_period;
        S_width = net_obj.S_width;
        zeta_thresh = net_obj.zeta_thresh;
        size_thresh = net_obj.size_thresh;
        force_thresh = net_obj.force_thresh;
        segma_thresh = net_obj.segma_thresh;
        fill_thresh = net_obj.fill_thresh;
        return *this;
    }

    //ceil function
    int Net::ceil(float N){
        return (N-int(N)==0)?int(N):int(N+1);
    }

    void Net::add_polygon(polygon *p)
    {
        polys.push_back(p);
    }
    void Net::remove_polygon(polygon *p)
    {
        typedef std::list<polygon*>::iterator poly_iter;
        poly_iter i;
        for(i = polys.begin(); i != polys.end(); ++i){
            if(*i == p)
                break;
        }
        delete (*i);
        polys.erase(i);
    }

    void Net::remove_vertex(vertex* v){
        typedef std::list<vertex*>::iterator ver_iter;
        deformable.remove(v);
        ver_iter i;
        for(i = verts.begin(); i != verts.end(); ++i){
            if(*i == v)
                break;
        }
        delete (*i);
        verts.erase(i);
    }

    //get neighbor of a polygon with edge v1-->v2
    //v1 precede v2 in the clockwise order
    //Prepared BY: Khaled M. Kamel
    polygon* Net::get_nbour(polygon* poly, vertex* v1, vertex* v2){
        typedef std::list<polygon*>::iterator poly_iter;
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;

        //for the case where there is no neighbor
        //the edge lies on the image border
        if( (v1->x==0&&v2->x==0) || (v1->x==image_width&&v2->x==image_width) ||
            (v1->y==0&&v2->y==0) || (v1->y==image_height&&v2->y==image_height) )
            return 0;

        //the fact that polygons are clockwise ordered in the vertex list of polygon
        //means that the neighbor comes after poly in one vertex
        //and before poly in the other vertex
        //find poly in both v1's and v2's polygon list
        poly_iter v1_poly, v2_poly, v1_next, v1_prev, v2_next, v2_prev;
        for(poly_iter i = v1->myPolygons.begin(); i != v1->myPolygons.end(); ++i){
            if((*i) == poly){
                v1_poly = i;
                break;
            }
        }
        for(poly_iter i = v2->myPolygons.begin(); i != v2->myPolygons.end(); ++i){
            if((*i) == poly){
                v2_poly = i;
                break;
            }
        }
        //get next and previous polygons in both v1 and v2
        v1_next = v1_poly;
        v1_next = (++v1_next == v1->myPolygons.end()) ? v1->myPolygons.begin() : v1_next;
        v1_prev = v1_poly;
        v1_prev = (--v1_prev == v1->myPolygons.end()) ? --v1->myPolygons.end() : v1_prev;
        v2_next = v2_poly;
        v2_next = (++v2_next == v2->myPolygons.end()) ? v2->myPolygons.begin() : v2_next;
        v2_prev = v2_poly;
        v2_prev = (--v2_prev == v2->myPolygons.end()) ? --v2->myPolygons.end() : v2_prev;

        //for the case where both v1 and v2 are only present in poly and neighbor
        if(*v1_next == *v2_prev && *v1_prev != *v2_next){
            return *v1_next;
        }
        else if(*v1_next != *v2_prev && *v1_prev == *v2_next){
            return *v2_next;
        }

        //for the case where both v1 and v2 are present in some other polygon
        //neighbor is the polygon that has v2---->v1 edge
        else if(*v1_next == *v2_prev && *v1_prev == *v2_next){
            //search for v2---->v1
            for(list_iter i = (*v1_next)->myVertices.begin(); i != (*v1_next)->myVertices.end(); ++i){
                for(ver_iter j = (*i).begin(); j != (*i).end(); ++j){
                    if(*j == v2){
                        //get next vertex
                        j = (++j == (*i).end()) ? (*i).begin() : j;
                        if(*j == v1)
                            return *v1_next;
                        else
                            return *v2_next;
                    }
                }
            }
        }
        return 0; //control should never reach this line
    }

    //chromaticity distance
    //by Shady R. Fam
    double Net::dist_chrom(color i,color j)
    {
        float theta=std::abs(i.H-j.H);
        if (theta>180) theta=360-theta;
        theta*=(3.1415/180);
        return sqrt(i.S*i.S + j.S*j.S - 2*i.S*j.S*cos(theta));
    }

    //calculate the color distance between a polygon average color and a pixel color
    //cylindrical distance
    //by Shady R. Fam
    double Net::dist_cyl(color i, color j)
    {

        float dI=(i.I-j.I);
        double dch=dist_chrom(i,j);
        return sqrt(dI*dI+dch*dch);
    }

    //Prepared BY: Mahmoud Kareem
    diplacement Net::measure_diplacement( vertex *v1,vertex *v2)
    {
        diplacement Result;

        if ( v2->x == v1->x ) //      Slope is not Determined ( value /0 )
            {
                     if ( v2->y > v1->y)
                         {
                             // deltaX= positive Value , deltaY= 0
                             Result.deltaY = 0.0;
                             Result.deltaX = S_width;  // displacement = + Sensitivity Region width

                         }
                     else
                         {
                                 // deltaX= Negative Value , deltaY= 0
                               Result.deltaY = 0.0;
                               Result.deltaX = -1* S_width;  // displacement = - Sensitivity Region width
                         }

            }
        else
            {
                // Edge Slope  = (  y2-y1 ) / ( x2-x1 )
             double Edge_slope=  double( v2->y - v1->y ) /double(  v2->x - v1->x );

                     if ( Edge_slope == 0 ) //straight Edge
                     {
                             if ( v2->x > v1->x)
                             {
                                 // deltaX=0 , deltaY= Negative value
                                 Result.deltaX = 0.0;
                                 Result.deltaY = -1*S_width;

                             }
                             else
                             {
                                   // deltaX=0 , deltaY= Positive value
                                   Result.deltaX = 0.0;
                                   Result.deltaY = S_width;
                             }

                     }
                     else if ( Edge_slope >0 ) // Positive Slope
                     {
                                 if ( v2->x > v1->x)
                                    {
                                        // deltaX= Postive value , deltaY= Negative value
                                        Result.deltaY = -1*sqrt(  ( double(S_width * S_width ) / ( (Edge_slope * Edge_slope)+1)) );
                                        Result.deltaX =  -1*Edge_slope* Result.deltaY;
                                    }
                             else
                                    {
                                      // deltaX= Negative  value , deltaY= Positive value
                                        Result.deltaY = sqrt(  ( double(S_width * S_width ) / ( (Edge_slope * Edge_slope)+1)) );
                                        Result.deltaX = -1*Edge_slope* Result.deltaY;
                                    }
                     }
                     else if ( Edge_slope <0 ) // Negative Slope
                     {
                                if ( v2->x > v1->x)
                                    {
                                    // deltaX=Negative  value , deltaY= Negative value
                                        Result.deltaY = -1*sqrt(  ( double(S_width * S_width ) / ( (Edge_slope * Edge_slope)+1)) );
                                        Result.deltaX =  -1*Edge_slope* Result.deltaY;

                                    }
                                else
                                    {
                                    // deltaX=Positive value , deltaY= Positive value
                                        Result.deltaY = sqrt(  ( double(S_width * S_width ) / ( (Edge_slope * Edge_slope)+1)) );
                                        Result.deltaX = -1* Edge_slope* Result.deltaY;
                                    }

                    }
            }
    return Result;
    }
//
//    // Return Senetivity Region Scanning  Dimentions  [ X_max , X_min  , Y_max  , Y_min ]
//    //Prepared BY: Mahmoud Kareem
//    scaning_d Net::measure_scaning_D ( vertex *v1 , vertex *v2 , vertex *v3 , vertex *v4)
//    {
//        scaning_d Result;
//        Result.X_max=v1->x;
//        Result.X_min=v1->x;
//        Result.Y_max=v1->y;
//        Result.Y_min=v1->y;
//        // search for max X && min X && Max Y , min Y
//        int temp_arrX[4]={ v1->x ,v2->x,v3->x,v4->x };
//        int temp_arrY[4]={ v1->y ,v2->y,v3->y,v4->y };
//
//        for ( int i= 1;i<4;i++)
//        {
//            // check for max and min X
//            if ( temp_arrX[i] > Result.X_max )
//                    Result.X_max =temp_arrX[i];
//            if ( temp_arrX[i] < Result.X_min )
//                    Result.X_min =temp_arrX[i];
//            // check for max and min Y
//            if ( temp_arrY[i] > Result.Y_max )
//                    Result.Y_max =temp_arrY[i];
//            if ( temp_arrY[i] < Result.Y_min )
//                    Result.Y_min =temp_arrY[i];
//        }
//
//        return Result;
//    }

    //check if The  given point ( x, y) is inside The Image Frame & if not Recordinate it to be inside The Frame
    //Prepared BY: Mahmoud Kareem
    void Net::coordinates_check( int * x, int *y, int max_X, int max_Y)
    {
       if ( *x >(max_X-1) )
            *x=(max_X-1);
       if ( *x <0 )
             *x=0;
       if ( *y >(max_Y-1) )
            *y=(max_Y-1);
       if ( *y <0)
            *y=0;
    }
//
//    //calculate the disparity for vertex [ V2 ] with both next vertex [ V3 ] and previous vertex [ V1 ]  in the polygon
//    //Prepared BY: Mahmoud Kareem
//     ver_disp Net::vertex_disparity(polygon* poly, vertex* v1, vertex* v2, vertex* v3,frame F)
//    {
//            /* initialize Results with 0 */
//        ver_disp result={0,0};
//
//         /***************************************************************************************************************************************
//         ****************************************** Calculate Disparity for [ V2 --------- V3 ] Edge ****************************************** */
//        /***************************************************************************************************************************************/
//        int x,y;
//
//         vertex center( int( ( v2->x + v3->x )/2 ) , int( ( v2->y + v3->y )/2 ) );   // Edge ( V2 ---------c---------V3)  Center Point
//
//        diplacement D= measure_diplacement( v2 , v3 );   // determine deltaX , deltaY
//
//         /* Vertex v after diplacement by deltaX , deltaY*/
//          x= int( v2->x +D.deltaX );
//          y= int( v2->y+D.deltaY );
//         coordinates_check( &x , &y ,F.width,F.height );      // check if points is inside The Image Frame
//          vertex v2_disp( x,y ) ;
//         /* Center Point after diplacement by deltaX , deltaY*/
//           x= int( center.x + D.deltaX );
//          y=  int( center.y+D.deltaY );
//         coordinates_check( &x , &y ,F.width,F.height );
//         vertex center_disp( x,y) ;
//
//         /* Sensitivity Region Polygon */
//
//         std::list<std::list<vertex*>> L;
//         std::list<vertex*> temp;
//         L.push_back(temp);
//         (*L.begin()).push_back(v2);
//         (*L.begin()).push_back(&v2_disp);
//         (*L.begin()).push_back(&center_disp);
//         (*L.begin()).push_back(&center);
//         polygon S_region(L);
//
////         /**debug**/
////         typedef std::list<std::list<vertex*>>::iterator list_iter;
////         typedef std::list<vertex*>::iterator ver_iter;
////         std::cout<<'\n'<<"polygon S_region : ";
////         for(list_iter m = S_region.myVertices.begin(); m != S_region.myVertices.end(); ++m){
////             for(ver_iter n = (*m).begin(); n != (*m).end(); ++n){
////                 std::cout<<(*n)->x<<","<<(*n)->y<<" ";
////             }
////             std::cout<<'\n';
////         }
////
////        std::ofstream file;
////        file.open("sensitivity.txt");
////        for(list_iter j = S_region.myVertices.begin(); j != S_region.myVertices.end(); ++j){
////            file << "(" <<int((*j).front()->x) << "," << int((*j).front()->y) << ") ";
////            for(ver_iter k = ++(*j).begin(); k != (*j).end(); ++k){
////                file << "(" << int((*k)->x) << "," << int((*k)->y) << ") ";
////                file << "(" << int((*k)->x) << "," << int((*k)->y) << ") ";
////            }
////            file << "(" <<int((*j).front()->x) << "," << int((*j).front()->y) << ")" << '\n';
////        }
////        file.close();
////
////         /**debug**/
//
//        /* Calculate Scanning Dimensions */
//         scaning_d DD=measure_scaning_D(v2 , &v2_disp , &center_disp , &center);
//
//         /* Counter to hold H value( number of pixel in the Sensetivity Region Which Achieve the Conditions ) */
//        int H=0;
//
////        /**debug**/
////        int count = 0;
////        /**debug**/
//
//        // Find Neighbour Polygon
//        polygon *nbour=get_nbour( poly , v2  ,  v3  );
//        //max disparity to hold the outer vertices on the border of the
//        if(nbour == 0)
//            result.H_next=0;
//        else{
//
////            /**debug**/
////            file.open("mine.txt");
////            int dist, nbourdist;
////            /**debug**/
//
//             // Scanning
//            for ( int y=DD.Y_min ; y<=DD.Y_max ; y++ )
//            {
//                for ( int x =DD.X_min ; x<=DD.X_max ; x++ )
//                {
//                    // check if The Pixel is Inside The Sensetivity Region Polygon and  in The Neghbour Polygon
//                    if ( nbour->PIP (  point(x,y) ) && S_region.PIP ( point(x,y)   )  )
//                    {
//
////                        /**debug**/
////                        file <<"("<< x <<","<< y <<") ";
////                        /**debug**/
////
////                        /**debug**/
////                        ++count;
////                        /**debug**/
//
//                                   // check if The Color Distance between Pixel Color and  its Polygon is > zeta_thresh
//                                   // and The Color Distance between Pixel Color and  The  Polygon Poly is <= zeta_thresh
//                                      // if This Conditions is True Do H++
//                              if (  dist_cyl( color( F.data[x][y] ) ,  nbour->c ) > zeta_thresh  && dist_cyl(color( F.data[x][y] )  , poly->c ) <= zeta_thresh ){
//                                                H++;
//
////                                                /**debug**/
////                                                nbourdist=dist_cyl( color( F.data[x][y] ) ,  nbour->c  );
////                                                dist=dist_cyl( color( F.data[x][y] ) ,  poly->c  );
////                                                file <<"("<< x <<","<< y <<") ";
////                                                /**debug**/
////                                                /**debug**/
////                                                std::ofstream colorfile;
////                                                colorfile.open("colorfile.txt");
////                                                colorfile << "nbour color H:"<<nbour->c.H<<" S:"<<nbour->c.S<<" I:"<<nbour->c.I<<'\n';
////                                                colorfile << "poly color H:"<<poly->c.H<<" S:"<<poly->c.S<<" I:"<<poly->c.I<<'\n';
////                                                colorfile << "pixel["<<x<<"]["<<y<<"]"<<" H:"<<color( F.data[x][y]).H<<" s:"<<color( F.data[x][y]).S<<" I:"<<color( F.data[x][y]).I<<'\n';
////                                                colorfile.close();
////                                                /**debug**/
//                              }
//
//                    }
//                }
//            }
//
////            /**debug**/
////            file.close();
////            std::cout<<'\n'<<"in: "<<count<<'\n';
////            /**debug**/
//
//            /* Store Result */
//            result.H_next=H;
//        }
//
//         /****************************************************************************************************************************************
//         ********************************************* Calculate Disparity for [  V1---------- V2 ] Edge************************************** */
//        /******************************************************************************************************************************************/
//
//          //vertex *v_prev=poly->prev_vertex(v);     /* get Previous Vertex*/
//
//          //( int( ( v->x + v_next->x )/2 ) , int( ( v->y + v_next->y )/2 ) );
//          // Edge ( Vprev---------c---------V )
//          center.x = int (  (v2->x + v1->x )/2 );
//          center.y = int (  (v2->y + v1->y )/2 ) ;
//
//         D= measure_diplacement( v1 , v2 );   // determine deltaX , deltaY
//
//         /* Vertex v after diplacement by deltaX , deltaY*/
//          x= int(v2->x +D.deltaX );
//          y=  int( v2->y +D.deltaY );
//        coordinates_check( &x , &y ,F.width,F.height );
//           v2_disp.x = x;
//           v2_disp.y = y;
//         /* Center Point after diplacement by deltaX , deltaY*/
//           x= int(center.x + D.deltaX );
//          y=  int( center.y+D.deltaY );
//        coordinates_check( &x , &y ,F.width,F.height );
//           center_disp.x=x;
//           center_disp.y=y;
//         /* Sensitivity Region Polygon */
//         std::list<std::list<vertex*>> L2;
//         std::list<vertex*> temp2;
//         L2.push_back(temp2);
//         (*L2.begin()).push_back(v2);
//         (*L2.begin()).push_back(&v2_disp);
//         (*L2.begin()).push_back(&center_disp);
//         (*L2.begin()).push_back(&center);
//         polygon S_region2(L2);
//
////         /**debug**/
////         std::cout<<'\n'<<"polygon S_region2 : ";
////         for(list_iter m = S_region2.myVertices.begin(); m != S_region2.myVertices.end(); ++m){
////             for(ver_iter n = (*m).begin(); n != (*m).end(); ++n){
////                 std::cout<<(*n)->x<<","<<(*n)->y<<" ";
////             }
////             std::cout<<'\n';
////         }
////
////         file.open("sensitivity.txt");
////         for(list_iter j = S_region2.myVertices.begin(); j != S_region2.myVertices.end(); ++j){
////             file << "(" <<int((*j).front()->x) << "," << int((*j).front()->y) << ") ";
////             for(ver_iter k = ++(*j).begin(); k != (*j).end(); ++k){
////                 file << "(" << int((*k)->x) << "," << int((*k)->y) << ") ";
////                 file << "(" << int((*k)->x) << "," << int((*k)->y) << ") ";
////             }
////             file << "(" <<int((*j).front()->x) << "," << int((*j).front()->y) << ")" << '\n';
////         }
////         file.close();
////
////         /**debug**/
//
//         /* Calculate Scanning Dimensions */
//          DD=measure_scaning_D(v2 ,&center,&center_disp, &v2_disp);
//
//         /* Counter to hold H value( number of pixel in the Sensetivity Region Which Achieve the Conditions ) */
//         H=0;
//
//
//
//         // Find Neighbour Polygon
//        nbour=get_nbour( poly ,v1 , v2);
//        if(nbour == 0)
//            result.H_prev=0;
//        else{
//
////            /**debug**/
////            file.open("mine.txt");
////            int dist, nbourdist;
////            count = 0;
////            /**debug**/
//
//            // Scanning
//            for ( int y=DD.Y_min ; y<=DD.Y_max ; y++ )
//            {
//                for ( int x =DD.X_min ; x<=DD.X_max ; x++ )
//                {
//                    // check if The Pixel is Inside The Sensetivity Region Polygon and  in The Neghbour Polygon
//                    if ( nbour->PIP ( point(x,y) ) && S_region2.PIP (  point(x,y)   )  )
//                    {
//
////                        /**debug**/
////                        file <<"("<< x <<","<< y <<") ";
////                        /**debug**/
////
////                         /**debug**/
////                        ++count;
////                        /**debug**/
//
//                               // check if The Color Distance between Pixel Color and  its Polygon is >zeta_thresh
//                             // and The Color Distance between Pixel Color and  The  Polygon Poly is <= zeta_thresh
//                           // if This Conditions is True Do H++
//                              if (  dist_cyl( color( F.data[x][y] ) ,  nbour->c  ) > zeta_thresh  && dist_cyl(color(F.data[x][y] )  ,poly->c ) <= zeta_thresh ){
//                                                H++;
//
////                                                /**debug**/
////                                                nbourdist=dist_cyl( color( F.data[x][y] ) ,  nbour->c  );
////                                                dist=dist_cyl( color( F.data[x][y] ) ,  poly->c  );
////                                                file <<"("<< x <<","<< y <<") ";
////                                                /**debug**/
////                                                /**debug**/
////                                                std::ofstream colorfile;
////                                                colorfile.open("colorfile.txt");
////                                                colorfile << "nbour color H:"<<nbour->c.H<<" S:"<<nbour->c.S<<" I:"<<nbour->c.I<<'\n';
////                                                colorfile << "poly color H:"<<poly->c.H<<" S:"<<poly->c.S<<" I:"<<poly->c.I<<'\n';
////                                                colorfile << "pixel["<<x<<"]["<<y<<"]"<<" H:"<<color( F.data[x][y]).H<<" s:"<<color( F.data[x][y]).S<<" I:"<<color( F.data[x][y]).I<<'\n';
////                                                colorfile.close();
////                                                /**debug**/
//                              }
//                    }
//                }
//            }
//
////            /**debug**/
////            file.close();
////            std::cout<<'\n'<<"in: "<<count<<'\n';
////            /**debug**/
//
//            /* Store Result */
//            result.H_prev=H;
//        }
//        return result;
//    }

    ver_disp Net::vertex_disparity(polygon* poly, vertex* v1, vertex* v2, vertex* v3, frame F){

        ver_disp result;
        point v2_point(v2->x, v2->y);
        //v1--->v2 mid point
        point v1v2_mid( (v1->x+v2->x)/2, (v1->y+v2->y)/2 );
        //v1--->v2 vector
        vector e1( v2->x-v1->x, v2->y-v1->y);
        //half v1--->v2 length
        int len = e1.mag()/2;
        //edge unit vector
        e1 = e1.unit_vec();
        //normal to edge unit vector
        vector v1v2_normal = e1.perpendicular();

        //initialize current pixel to non-existent one
        point pixel;
        result.H_prev = 0;
        polygon* nbour = get_nbour(poly, v1, v2);
        //compute H_prev
        for(int j = 1; j <= S_width; ++j){
            for(int i = 0; i < len; ++i){
                pixel = v1v2_mid.move(v1v2_normal*j+ e1*i);
                if(pixel.x>=0 && pixel.x<image_width && pixel.y>=0 && pixel.y<image_height &&
                   //nbour->PIP(pixel) &&
                   //color distance between current pixel and neighbor above color theshold
                   (dist_cyl( F.data[pixel.x][pixel.y], nbour->c) > zeta_thresh) &&
                   //color distance between current pixel and poly below or equal color theshold
                   (dist_cyl( F.data[pixel.x][pixel.y], poly->c) <= zeta_thresh) ){

                    ++result.H_prev;
                }
            }
        }

        //v2--->v3 vector
        vector e2( v3->x-v2->x, v3->y-v2->y);
        //half v2--->v3 length
        len = e2.mag()/2;
        //edge unit vector
        e2 = e2.unit_vec();
        //normal to edge unit vector
        vector v2v3_normal = e2.perpendicular();

        result.H_next = 0;
        nbour = get_nbour(poly, v2, v3);
        //compute H_next
        for(int j = 1; j <= S_width; ++j){
            for(int i = 1; i <= len; ++i){
                pixel = v2_point.move(v2v3_normal*j+ e2*i);
                if(pixel.x>=0 && pixel.x<image_width && pixel.y>=0 && pixel.y<image_height &&
                   //nbour->PIP(pixel) &&
                   //color distance between current pixel and neighbor above color theshold
                   (dist_cyl( F.data[pixel.x][pixel.y], nbour->c) > zeta_thresh) &&
                   //color distance between current pixel and poly below or equal color theshold
                   (dist_cyl( F.data[pixel.x][pixel.y], poly->c) <= zeta_thresh) ){

                    ++result.H_next;
                }
            }
        }

        return result;
    }

    //go through all vertices in all polygons
    //calculate the deformation forces and update the location of the vertices
    //completing a full deformation cycle
    //return the maximum vertex displacement in one full cycle
    //Prepared BY : Khaled M. Kamel
    float Net::deformation(frame F){
        typedef std::list<polygon*>::iterator poly_iter;
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;
        vertex *current , *next , *prev;
        //disparity for a vertex in a polygon
        bool found;
        ver_disp D;
        //edge length
        float len;
        //force magnitude
        float mag;
        //perpendicular to edge unit vectors
        //perpendicular force vector
        vector per_unit, force_vec;
        //maximum vertex displacement
        float max_disp = 0;
        //go through all vertices in all polygons
        for(poly_iter i = polys.begin(); i != polys.end(); ++i){
            //for each vertex k in polygon i
            for(list_iter j = (*i)->myVertices.begin(); j != (*i)->myVertices.end(); ++j){
                for(ver_iter k = (*j).begin(); k != (*j).end(); ++k){
                    current = *k;
                    found = false;
                    for(ver_iter m = deformable.begin(); m != deformable.end(); ++m){
                        if(current == *m){
                            found = true;
                            break;
                        }
                    }
                    if(found){
                        /*for all vertices that are not on the border of the image that have moved
                          in the previous deformation cycle*/
                        if( current->x!=0 && current->x!=image_width && current->y!=0 && current->y!=image_height ){
                            //if current vertex is the last one then next vertex is the first one (circular list)
                            next = (++k == (*j).end()) ? *(*j).begin() : *(k);
                            --k;
                            //if current vertex is the first one then the previous vertex is the last one (circular list)
                            prev = (--k == (*j).end()) ? *(--(*j).end()) : *(k);
                            ++k;
                            //calculate the disparity for the current vertex with both next and previous vertices
                            D = vertex_disparity(*i, prev, current, next, F);
                            if( (current->x != next->x) || (current->y != next->y) ){
                                //calculate the force on the current vertex in edge with the next vertex
                                //calculate the edge length
                                len = sqrt( (next->x - current->x)*(next->x - current->x) +
                                           (next->y - current->y)*(next->y - current->y) );
                                //calculate the force magnitude
                                mag = D.H_next / len;
                                //a unit vector perpendicular to the edge (current ---> next)
                                per_unit.update_coords((next->y - current->y)/len, -(next->x - current->x)/len);
                                //the force vector
                                force_vec = per_unit * mag;
                                //update total force on vertex
                                current->total_force += force_vec;
                            }
                            //---------------------------------------------
                            if( (current->x != prev->x) || (current->y != prev->y) ){
                                //calculate the force on the current in edge with the previous vertex
                                //calculate the edge length
                                len = sqrt( (prev->x - current->x)*(prev->x - current->x) +
                                           (prev->y - current->y)*(prev->y - current->y) );
                                //calculate the force magnitude
                                mag = D.H_prev / len;
                                //a unit vector perpendicular to the edge (prev ---> current)
                                per_unit.update_coords((current->y - prev->y)/len, -(current->x - prev->x)/len);
                                //the force vector
                                force_vec = per_unit * mag;
                                //update total force on vertex
                                current->total_force += force_vec;
                            }
                        }
                    }
                }
            }
        }
        //update all net vertices to the new locations
        vector disp_vec;
        int ind, size;
        bool intersect;
        ver_iter v_next;
        float a1,b1,c1,a2,b2,c2;
        float s, t;
        float denom;
        for(std::list<vertex*>::iterator i = verts.begin(); i != verts.end(); ++i){
            disp_vec.update_coords((*i)->x, (*i)->y);
            //move current vertex to the new location
            (*i)->update_loc( (*i)->x+(*i)->total_force.x, (*i)->y+(*i)->total_force.y );
            if((*i)->x == disp_vec.x && (*i)->y == disp_vec.y){
                deformable.remove(*i);
            }
            else{
                size = (*i)->myPolygons.size();
                ind = 0;
                vertex* nbour_verts[size];
                for(poly_iter j = (*i)->myPolygons.begin(); j != (*i)->myPolygons.end(); ++j){
                    for(list_iter k = (*j)->myVertices.begin(); k != (*j)->myVertices.end(); ++k){
                        ver_iter v;
                        for(v = (*k).begin(); v != (*k).end(); ++v){
                            if(*v == *i){
                                v = (--v == (*k).end()) ? --(*k).end() : v;
                                nbour_verts[ind] = *v;
                                ++ind;
                                break;
                            }
                        }
                        if(v != (*k).end())
                            break;
                    }
                }
                intersect = false;
                for(poly_iter j = (*i)->myPolygons.begin(); j != (*i)->myPolygons.end(); ++j){
                    for(list_iter k = (*j)->myVertices.begin(); k != (*j)->myVertices.end(); ++k){
                        for(ver_iter v = (*k).begin(); v != (*k).end(); ++v){
                            v_next = v;
                            v_next = (++v_next == (*k).end()) ? (*k).begin() : v_next;
                            if(*v != *i && *v_next != *i){
                                for(ind = 0; ind < size; ++ind){
                                    if(*v != nbour_verts[ind] && *v_next != nbour_verts[ind]){
                                        a1 = nbour_verts[ind]->x - (*i)->x;
                                        b1 = (*v)->x - (*v_next)->x;
                                        c1 = (*v)->x - (*i)->x;
                                        a2 = nbour_verts[ind]->y - (*i)->y;
                                        b2 = (*v)->y - (*v_next)->y;
                                        c2 = (*v)->y - (*i)->y;
                                        denom = a1*b2 - a2*b1;
                                        if(denom != 0){
                                            s = (c1*b2 - c2*b1) / denom;
                                            t = (a1*c2 - a2*c1) / denom;
                                            if(s>=0 && s<=1 && t>=0 && t<=1 ){
                                                intersect = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                if(intersect)
                                    break;
                            }
                        }
                        if(intersect)
                            break;
                    }
                    if(intersect)
                            break;
                }
                if(intersect){
                    (*i)->update_loc(disp_vec.x, disp_vec.y);
                }
                else{
                    disp_vec.update_coords( (*i)->x-disp_vec.x, (*i)->y-disp_vec.y );
                    max_disp = std::max(max_disp, disp_vec.mag());
                }
            }
            //reset the total force on current vertex
            (*i)->total_force.update_coords(0,0);

        }
        return max_disp;
    }

    polygon* Net::merge_utility(polygon* p1, polygon* p2, frame F){

        typedef std::list<polygon*>::iterator poly_iter;
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;

        std::list<std::list<vertex*>> merged_poly;
        polygon* p3 = new polygon(merged_poly);
        //add the merged polygon to the net
        polys.push_back(p3);
        /*  shared_verts is a list of shared connected_verts
            connected_verts is a list of vertices iterators{V1 , V2 , ... , Vn} where :
                Vi is a shared vertex
                Vi->Vi+1 is an edge in p1
                Vi+1->Vi is a edge in p2
                V0->V1 is an edge in p1
                V1->V0 is not an edge in p2
        */
        typedef std::list<ver_iter> connected_verts;
        typedef std::list<connected_verts> shared_verts;
        shared_verts p1_sv;
        shared_verts p2_sv;
        connected_verts cv;

        list_iter p1_outerframe;
        list_iter p2_sharedlist;
        //save iterators to shared vertices in both p1 and p2
        for(int i = 0; i < 2; ++i){
            p1_outerframe = p1->myVertices.begin();
            p2_sharedlist = p2->myVertices.end();
            ver_iter p1_current;
            ver_iter p2_current;
            //find the shared list in p2
            for(p1_current = (*p1_outerframe).begin(); p1_current != (*p1_outerframe).end(); ++p1_current){
                for(list_iter j = p2->myVertices.begin(); j != p2->myVertices.end(); ++j){
                    for(p2_current = (*j).begin(); p2_current != (*j).end(); ++p2_current){
                        if(*p2_current == *p1_current){
                            p2_sharedlist = j;
                            break;
                        }
                    }
                    if(p2_sharedlist != p2->myVertices.end())
                        break;
                }
                if(p2_sharedlist != p2->myVertices.end())
                    break;
            }
            //switch p1 and p2 if no shared vertices are found in p1's outer frame
            if(p2_sharedlist == p2->myVertices.end()){
                polygon* temp = p1;
                p1 = p2;
                p2 = temp;
                continue;
            }

            //find V1
            ver_iter p1_prev, p2_next, temp;
            while(p1_current != (*p1_outerframe).end()){
                //find current vertex in p2
                for(p2_current = (*p2_sharedlist).begin(); p2_current != (*p2_sharedlist).end(); ++p2_current){
                    if(*p2_current == *p1_current)
                        break;
                }
                //go to the next vertex in p1 if current vertex is not shared
                if(p2_current == (*p2_sharedlist).end()){
                    ++p1_current;
                    continue;
                }
                temp = p1_current;
                p1_prev = (--temp == (*p1_outerframe).end()) ? --(*p1_outerframe).end() : temp;
                temp = p2_current;
                p2_next = (++temp == (*p2_sharedlist).end()) ? (*p2_sharedlist).begin() : temp;
                if(*p1_prev != *p2_next){
                    //V1 found
                    break;
                }
                ++p1_current;
            };

            //if V1 is not found
            //all edges of p1 are shared with p2 (p1 is a hole in p2)
            if(p1_current == (*p1_outerframe).end()){
                //form the merged polygon
                //the new outer frame is p2's outer frame
                //the new holes are p2 holes minus p1 hole in p2 and all p1 holes
                p2->myVertices.erase(p2_sharedlist);
                p3->myVertices = p2->myVertices;
                p3->myVertices.insert(p3->myVertices.end(), ++p1->myVertices.begin(), p1->myVertices.end());

                //clean up
                //remove p1's outer frame from the net
                for(ver_iter j = (*p1_outerframe).begin(); j != (*p1_outerframe).end(); ++j){
                    remove_vertex(*j);
                }
                //for all vertices of p1 holes
                //remove p1 and add p3 at exact same location in list
                for(list_iter j = ++(p1->myVertices.begin()); j != p1->myVertices.end(); ++j){
                    for(ver_iter k = (*j).begin(); k != (*j).end(); ++k){
                        for(poly_iter m = (*k)->myPolygons.begin(); m != (*k)->myPolygons.end(); ++m){
                            if(*m == p1){
                                (*k)->myPolygons.insert( (*k)->myPolygons.erase(m), p3 );
                                break;
                            }
                        }
                    }
                }
                //for all vertices of p2 (except p1 hole that has already been removed)
                //remove p2 and add p3 at exact same location in list
                for(list_iter j = p2->myVertices.begin(); j != p2->myVertices.end(); ++j){
                    for(ver_iter k = (*j).begin(); k != (*j).end(); ++k){
                        for(poly_iter m = (*k)->myPolygons.begin(); m != (*k)->myPolygons.end(); ++m){
                            if(*m == p2){
                                (*k)->myPolygons.insert( (*k)->myPolygons.erase(m), p3 );
                                break;
                            }
                        }
                    }
                }
                p3->c.R = (p1->c.R * p1->area + p2->c.R * p2->area)/(p1->area + p2->area);
                p3->c.G = (p1->c.G * p1->area + p2->c.G * p2->area)/(p1->area + p2->area);
                p3->c.B = (p1->c.B * p1->area + p2->c.B * p2->area)/(p1->area + p2->area);
                computeHSI(p3->c);
                p3->area = p1->area + p2->area;
                //remove p1 and p2 from the net
                remove_polygon(p1);
                remove_polygon(p2);
                return p3;
            }

            //form the shared lists of p1 and p2
            p1_sv.push_back(cv);
            shared_verts::iterator p1_currset = p1_sv.begin();
            (*p1_currset).push_back(p1_current);
            p2_sv.push_back(cv);
            shared_verts::iterator p2_currset = p2_sv.begin();
            (*p2_currset).push_back(p2_current);
            ver_iter p1_start = p1_current;
            temp = p1_current;
            p1_current = (++temp == (*p1_outerframe).end()) ? (*p1_outerframe).end() : temp;
            while(p1_current != p1_start){
                //find current vertex in p2
                for(p2_current = (*p2_sharedlist).begin(); p2_current != (*p2_sharedlist).end(); ++p2_current){
                    if(*p2_current == *p1_current)
                        break;
                }
                //go to the next vertex in p1 if current vertex is not shared
                if(p2_current == (*p2_sharedlist).end()){
                    p1_current = (++p1_current == (*p1_outerframe).end()) ? (*p1_outerframe).end() : p1_current;
                    continue;
                }
                temp = p1_current;
                p1_prev = (--temp == (*p1_outerframe).end()) ? --(*p1_outerframe).end() : temp;
                temp = p2_current;
                p2_next = (++temp == (*p2_sharedlist).end()) ? (*p2_sharedlist).begin() : temp;
                if(*p1_prev != *p2_next){
                    p1_sv.push_back(cv);
                    p1_currset = --p1_sv.end();
                    p2_sv.push_back(cv);
                    p2_currset = --p2_sv.end();
                }
                (*p1_currset).push_back(p1_current);
                (*p2_currset).push_back(p2_current);
                p1_current = (++p1_current == (*p1_outerframe).end()) ? (*p1_outerframe).end() : p1_current;
            }
            break;
        }

        //form subpolygons that form the merged polygon
        /*  for shared_verts {connected_verts 1 , connected_verts 2, ... , connected_verts n}
            where connected_verts i {V1, V2 , ... , Vn}
            a subpolygon is formed in the following fashion:
                subpoly i = vertices of p1 from connected_verts i::Vn to connected_verts i+1::V1
                          + vertices of p2 from connected_verts i+1::V1 to connected_verts i::Vn
            treating the shared_verts as a cylindrical list (connected_verts n+1 = connected_verts 1)
        */
        std::list<vertex*> subpolys[p1_sv.size()];
        shared_verts::iterator temp;
        shared_verts::iterator p1_currset = p1_sv.begin();
        temp = p1_currset;
        shared_verts::iterator p1_nextset = (++temp == p1_sv.end()) ? p1_sv.begin() : temp;
        shared_verts::iterator p2_currset = p2_sv.begin();
        temp = p2_currset;
        shared_verts::iterator p2_nextset = (++temp == p2_sv.end()) ? p2_sv.begin() : temp;
        ver_iter currvert;
        int p1_sv_size = p1_sv.size();
        for(int i = 0; i < p1_sv_size; ++i){
            currvert = (*p1_currset).back();
            while(currvert != (*p1_nextset).front()){
                subpolys[i].push_back(*currvert);
                currvert = (++currvert == (*p1_outerframe).end()) ? (*p1_outerframe).begin() : currvert;
            }
            currvert = (*p2_nextset).front();
            while(currvert != (*p2_currset).back()){
                subpolys[i].push_back(*currvert);
                currvert = (++currvert == (*p2_sharedlist).end()) ? (*p2_sharedlist).begin() : currvert;
            }
            p1_currset = p1_nextset;
            p1_nextset = (++p1_nextset == p1_sv.end()) ? p1_sv.begin() : p1_nextset;
            p2_currset = p2_nextset;
            p2_nextset = (++p2_nextset == p2_sv.end()) ? p2_sv.begin() : p2_nextset;
        }

        /*  two subpolygons might end up sharing a single vertex.
            this can be easily detected when a connected_verts list contains only a single vertex iterator.
            the shared vertex should be duplicated.
        */
        ver_iter p1_nextver , p1_prevver;
        polygon* p1_nbour;
        int currsubpoly = 0;
        int prevsubpoly = p1_sv.size()-1;
        for(p1_currset = p1_sv.begin(); p1_currset != p1_sv.end(); ++p1_currset){
            //if a connected_verts list contains only a single vertex iterator
            if((*p1_currset).size() == 1){
                currvert = (*p1_currset).front();
                temp = p1_currset;
                p1_nextset = (++temp == p1_sv.end()) ? p1_sv.begin() : temp;
                //find a neighbor of p1 that is part of subpolys[currindex]
                p1_nextver = currvert;
                p1_nextver = (++p1_nextver == (*p1_outerframe).end()) ? (*p1_outerframe).begin() : p1_nextver;
                p1_nbour = get_nbour(p1, *currvert, *p1_nextver);
                //split the list of neighbors in the shared vertex into the two groups enclosed between p1 and p2
                //identifying the groups that contains p1_nbour in the process
                std::list<polygon*> groupone;
                std::list<polygon*> grouptwo;
                //true if group one contains p1_nbour
                bool nbourinone = true;
                //find p1 or p2
                poly_iter currpoly;
                for(currpoly = (*currvert)->myPolygons.begin(); currpoly != (*currvert)->myPolygons.end(); ++currpoly){
                    if( (*currpoly == p1) || (*currpoly == p2) )
                        break;
                }
                //form the two groups skipping p1 and p2
                currpoly = (++currpoly == (*currvert)->myPolygons.end()) ? (*currvert)->myPolygons.begin() : currpoly;
                while( (*currpoly != p1) && (*currpoly != p2)){
                    groupone.push_back(*currpoly);
                    currpoly = (++currpoly == (*currvert)->myPolygons.end()) ? (*currvert)->myPolygons.begin() : currpoly;
                }
                currpoly = (++currpoly == (*currvert)->myPolygons.end()) ? (*currvert)->myPolygons.begin() : currpoly;
                while( (*currpoly != p1) && (*currpoly != p2)){
                    grouptwo.push_back(*currpoly);
                    if(*currpoly == p1_nbour)
                        nbourinone = false;
                    currpoly = (++currpoly == (*currvert)->myPolygons.end()) ? (*currvert)->myPolygons.begin() : currpoly;
                }
                //add p3 to both groups
                groupone.push_back(p3);
                grouptwo.push_back(p3);
                //identify the group that have fewer and larger number of polygons
                std::list<polygon*>* fewer = &groupone;
                std::list<polygon*>* larger = &grouptwo;
                bool onehasfewer = true;
                if(groupone.size() < grouptwo.size()){
                    fewer = &grouptwo;
                    larger = &groupone;
                    onehasfewer = false;
                }
                //create a new copy of the shared vertex with the group the has fewer number of polygons
                vertex* currdouble = new vertex((*currvert)->x, (*currvert)->y, *fewer);
                verts.push_back(currdouble);
                deformable.push_back(currdouble);
                //update the polygon list of the shared vertex with the group that has larger number of polygons
                (*currvert)->myPolygons = *larger;
                //remove the shared vertex from the the polygon of fewer group and the corresponding subpoly
                //and add the new copy at the exact same location in the vertex list of each polygon
                list_iter outerframe;
                for(currpoly = fewer->begin(); currpoly != --(fewer->end()); ++currpoly){
                    outerframe = (*currpoly)->myVertices.begin();
                    for(ver_iter v = (*outerframe).begin(); v != (*outerframe).end(); ++v){
                        if(*v == *currvert){
                            (*outerframe).insert((*outerframe).erase(v), currdouble);
                            break;
                        }
                    }
                }
                if(nbourinone == onehasfewer){
                    for(ver_iter v = subpolys[currsubpoly].begin(); v != subpolys[currsubpoly].end(); ++v){
                        if(*v == *currvert){
                            subpolys[currsubpoly].insert(subpolys[currsubpoly].erase(v), currdouble);
                            break;
                        }
                    }
                }
                else{
                    for(ver_iter v = subpolys[prevsubpoly].begin(); v != subpolys[prevsubpoly].end(); ++v){
                        if(*v == *currvert){
                            subpolys[prevsubpoly].insert(subpolys[prevsubpoly].erase(v), currdouble);
                            break;
                        }
                    }
                }

            }
            prevsubpoly = currsubpoly;
            ++currsubpoly;
        }

        //form the merged polygon
        //if shared edges are not in p2's outer frame
        if(p2_sharedlist != p2->myVertices.begin()){
            //the new outer frame is p2's outer frame
            //the new holes are p2 holes minus p1 hole in p2 , all p1 holes and all subpolys
            p3->myVertices = p2->myVertices;
            p3->myVertices.remove(*p2_sharedlist);
            p3->myVertices.insert(p3->myVertices.end(), ++p1->myVertices.begin(), p1->myVertices.end());
            for(int i = 0; i < p1_sv_size; ++i){
                p3->myVertices.push_back(subpolys[i]);
            }
        }
        //if shared edges are in p2's outer frame
        else{
            //the new outer frame is the subpoly that clockwise orinted
            //the new holes are all other subpolys , all p1 holes and all p2 holes
            //find the new outer frame
            int ind = 0;
            //determinant of orientation matrix
            float orientDet = -1;
            while(orientDet < 0){
                //find the vertices that lie on the convex hull of the current subpoly
                int minx = F.width, maxx = 0, miny = F.height, maxy = 0;
                for(ver_iter j = subpolys[ind].begin(); j != subpolys[ind].end(); ++j){
                    minx = ((*j)->x < minx) ? (*j)->x :minx;
                    maxx = ((*j)->x > maxx) ? (*j)->x :maxx;
                    miny = ((*j)->y < miny) ? (*j)->y :miny;
                    maxy = ((*j)->y > maxy) ? (*j)->y :maxy;
                }
                std::list<ver_iter> smallestx, largestx, smallesty, largesty;
                for(ver_iter j = subpolys[ind].begin(); j != subpolys[ind].end(); ++j){
                    if((*j)->x == minx){
                        smallestx.push_back(j);
                        continue;
                    }
                    if((*j)->x == maxx){
                        largestx.push_back(j);
                        continue;
                    }
                    if((*j)->y == miny){
                        smallesty.push_back(j);
                        continue;
                    }
                    if((*j)->y == maxy){
                        largesty.push_back(j);
                    }
                }
                std::list<ver_iter> onhull;
                miny = F.height;
                maxy = 0;
                for(std::list<ver_iter>::iterator j = smallestx.begin(); j != smallestx.end(); ++j){
                    miny = ((*(*j))->y < miny) ? (*(*j))->y :miny;
                    maxy = ((*(*j))->y > maxy) ? (*(*j))->y :maxy;
                }
                for(std::list<ver_iter>::iterator j = smallestx.begin(); j != smallestx.end(); ++j){
                    if( ((*(*j))->y == miny) || ((*(*j))->y == maxy) )
                        onhull.push_back(*j);
                }
                miny = F.height;
                maxy = 0;
                for(std::list<ver_iter>::iterator j = largestx.begin(); j != largestx.end(); ++j){
                    miny = ((*(*j))->y < miny) ? (*(*j))->y :miny;
                    maxy = ((*(*j))->y > maxy) ? (*(*j))->y :maxy;
                }
                for(std::list<ver_iter>::iterator j = largestx.begin(); j != largestx.end(); ++j){
                    if( ((*(*j))->y == miny) || ((*(*j))->y == maxy) )
                        onhull.push_back(*j);
                }
                minx = F.width;
                maxx = 0;
                for(std::list<ver_iter>::iterator j = smallesty.begin(); j != smallesty.end(); ++j){
                    minx = ((*(*j))->x < minx) ? (*(*j))->x :minx;
                    maxx = ((*(*j))->x > maxx) ? (*(*j))->x :maxx;
                }
                for(std::list<ver_iter>::iterator j = smallesty.begin(); j != smallesty.end(); ++j){
                    if( ((*(*j))->x == minx) || ((*(*j))->x == maxx) )
                        onhull.push_back(*j);
                }
                minx = F.width;
                maxx = 0;
                for(std::list<ver_iter>::iterator j = largesty.begin(); j != largesty.end(); ++j){
                    minx = ((*(*j))->x < minx) ? (*(*j))->x :minx;
                    maxx = ((*(*j))->x > maxx) ? (*(*j))->x :maxx;
                }
                for(std::list<ver_iter>::iterator j = largesty.begin(); j != largesty.end(); ++j){
                    if( ((*(*j))->x == minx) || ((*(*j))->x == maxx) )
                        onhull.push_back(*j);
                }
                std::list<ver_iter>::iterator currentvert = onhull.begin();
                ver_iter prevvert, nextvert;
                orientDet = 0;
                while(orientDet ==0){
                    prevvert = *currentvert;
                    prevvert = (--prevvert == subpolys[ind].end()) ? --subpolys[ind].end() : prevvert;
                    nextvert = *currentvert;
                    nextvert = (++nextvert == subpolys[ind].end()) ? subpolys[ind].begin() : nextvert;
                    orientDet = (*(*currentvert))->x * (*nextvert)->y
                              - (*(*currentvert))->y * (*nextvert)->x
                              + (*prevvert)->x * (*(*currentvert))->y
                              - (*prevvert)->y * (*(*currentvert))->x
                              + (*prevvert)->y * (*nextvert)->x
                              - (*prevvert)->x * (*nextvert)->y;
                    ++currentvert;
                }
                ++ind;
            }
            int p3outer_frame = ind-1;
            //add the new outer frame
            p3->myVertices.push_back(subpolys[p3outer_frame]);
            //add all other subpolys
            for(int i = 0; i < p1_sv_size; ++i){
                if(i != p3outer_frame)
                   p3->myVertices.push_back(subpolys[i]);
            }
            //add all p1's holes
            p3->myVertices.insert(p3->myVertices.end(), (++p1->myVertices.begin()), p1->myVertices.end());
            //add all p2's holes
            p3->myVertices.insert(p3->myVertices.end(), (++p2->myVertices.begin()), p2->myVertices.end());
        }

        /**cleanup**/
        /*
            divide p1 and p2 vertices into 4 categories
            1 - non shared vertices of p1 (remove p1 and add p3 at exact same location in list)
            2 - non shared vertices of p2 (remove p2 and add p3 at exact same location in list)
            3 - endpoints of connected sets of shared vertices (remove p1 and p2 and add p3 at exact same location in list)
            4 - non endpoints of connected sets of shared vertices (removed totally from the net)
        */
        //for category 4 : remove totally from the net
        for(shared_verts::iterator i = p1_sv.begin(); i != p1_sv.end(); ++i){
            if((*i).size() > 2){
                for(connected_verts::iterator j = ++(*i).begin(); j != --(*i).end(); ++j){
                     remove_vertex(*(*j));
                }
            }
        }
        //for category 3 : remove p1 and p2 and add p3 at exact same location in list
        vertex *first, *last;
        for(shared_verts::iterator i = p1_sv.begin(); i != p1_sv.end(); ++i){
            if((*i).size() != 1){
                first = *(*i).front();
                first->myPolygons.remove(p1);
                for(poly_iter k = first->myPolygons.begin(); k != first->myPolygons.end(); ++k){
                    if(*k == p2){
                        first->myPolygons.insert( first->myPolygons.erase(k) , p3 );
                        break;
                    }
                }
                last = *(*i).back();
                last->myPolygons.remove(p1);
                for(poly_iter k = last->myPolygons.begin(); k != last->myPolygons.end(); ++k){
                    if(*k == p2){
                        last->myPolygons.insert( last->myPolygons.erase(k) , p3 );
                        break;
                    }
                }
            }
        }
        //remove shared vertices from p1 , p2
        //to form categories 1 and 2
        for(shared_verts::iterator i = p1_sv.begin(); i != p1_sv.end(); ++i){
            for(connected_verts::iterator j = (*i).begin(); j != (*i).end(); ++j){
                (*p1_outerframe).erase(*j);
            }
        }
        for(shared_verts::iterator i = p2_sv.begin(); i != p2_sv.end(); ++i){
            for(connected_verts::iterator j = (*i).begin(); j != (*i).end(); ++j){
                (*p2_sharedlist).erase(*j);
            }
        }
        //for category 1 : remove p1 and add p3 at exact same location in list
        for(list_iter i = p1->myVertices.begin(); i != p1->myVertices.end(); ++i){
            for(ver_iter j = (*i).begin(); j != (*i).end(); ++j){
                for(poly_iter k = (*j)->myPolygons.begin(); k != (*j)->myPolygons.end(); ++k){
                    if(*k == p1){
                        (*j)->myPolygons.insert( (*j)->myPolygons.erase(k) , p3 );
                        break;
                    }
                }
            }
        }
        //for category 2 : remove p2 and add p3 at exact same location in list
        for(list_iter i = p2->myVertices.begin(); i != p2->myVertices.end(); ++i){
            for(ver_iter j = (*i).begin(); j != (*i).end(); ++j){
                for(poly_iter k = (*j)->myPolygons.begin(); k != (*j)->myPolygons.end(); ++k){
                    if(*k == p2){
                        (*j)->myPolygons.insert( (*j)->myPolygons.erase(k) , p3 );
                        break;
                    }
                }
            }
        }
        p3->c.R = (p1->c.R * p1->area + p2->c.R * p2->area)/(p1->area + p2->area);
        p3->c.G = (p1->c.G * p1->area + p2->c.G * p2->area)/(p1->area + p2->area);
        p3->c.B = (p1->c.B * p1->area + p2->c.B * p2->area)/(p1->area + p2->area);
        computeHSI(p3->c);
        p3->area = p1->area + p2->area;
        //remove p1 and p2 from the net
        remove_polygon(p1);
        remove_polygon(p2);
        return p3;
    }


    //merge neighboring polygons having average color difference below 60
    //Prepared BY: Khaled M. Kamel
    void Net::color_merge(frame F){
        typedef std::list<polygon*>::iterator poly_iter;
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;

        //create a list of pointer to all polygons(to iterate through)
        std::list<polygon*> temp_polys;
        for(poly_iter i = polys.begin(); i != polys.end(); ++i)
            temp_polys.push_back(*i);

        polygon* i;
        list_iter first_list;
        vertex *v1, *v2;
        polygon *nbour, *nearest_nbour;
        double color_dist, nearest_color_dist;

        //go through all polygons
        while(!temp_polys.empty()){
            i = *temp_polys.begin();
            //initialize color distance of nearest neighbor
            nearest_color_dist = 100;
            //find the nearest neighbor in color
            //go through all outer edges of polygon i
            first_list = i->myVertices.begin();
            for(ver_iter k = (*first_list).begin(); k != (*first_list).end(); ++k){
                v1 = *k;
                v2 = (++k == (*first_list).end()) ? (*first_list).front() : *k;
                --k;
                //get the neighbor of polygon i with edge v1--->v2
                nbour = get_nbour(i, v1, v2);
                //skip the edge if it has no neighbor
                if(nbour == 0)
                    continue;
                //calculate the color distance of this neighbor
                color_dist = dist_cyl(i->c,nbour->c);
                //update the nearest neighbor in color
                if (color_dist <= 60 && color_dist < nearest_color_dist){
                    nearest_color_dist = color_dist;
                    nearest_nbour = nbour;
                }
            }
            if (nearest_color_dist <= zeta_thresh){
                //remove neighbor and add the merged polygon
                temp_polys.push_back(merge_utility(i, nearest_nbour, F));
                temp_polys.remove(nearest_nbour);
            }
            //remove polygon i
            temp_polys.pop_front();
        }
    }


    //merge a polygon with an area <= size threshold to a neighboring polygon having the nearest average color
    //Prepared BY: Khaled M. Kamel
    void Net::size_merge(frame F){
        typedef std::list<polygon*>::iterator poly_iter;
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;

        /*make a list of all polygons with area <= size threshold*/
        std::list<polygon*> small_polys;
        //go through all polygons
        for(poly_iter i = polys.begin(); i != polys.end(); ++i){
            //save pointer to polygons with area <= size threshold
            if((*i)->area <= size_thresh)
                small_polys.push_back(*i);
        }

        vertex *v1, *v2;
        list_iter first_list;
        //nearest neighbor in color
        polygon* nearest_nbour;
        polygon* nbour;
        //color distance of nearest neighbor
        double nearest_color_dist;
        double color_dist;
        polygon* i;
        /*merge each polygon with its neighboring polygon having the nearest average color*/
        while(!small_polys.empty()){
            i = *small_polys.begin();
            //check the size for the newly added merged polygons, as they might have an area <= size threshold
            if(i->area <= size_thresh){
                //initialize nearest neighbor
                //find an edge with a neighbor
                bool found = false;
                first_list = i->myVertices.begin();
                for(ver_iter k = (*first_list).begin(); k != (*first_list).end(); ++k){
                    v1 = *k;
                    v2 = (++k == (*first_list).end()) ? (*first_list).front() : *k;
                    --k;
                    nearest_nbour = get_nbour(i, v1, v2);
                    //break if found a neighbor
                    if(nearest_nbour != 0){
                        found = true;
                        break;
                    }
                }
                if(!found){
                    small_polys.pop_front();
                    continue;
                }
                //initialize color distance of nearest neighbor
                nearest_color_dist = dist_cyl(i->c, nearest_nbour->c);
                //get the nearest neighbor in color
                //go through all edges of polygon i
                for(ver_iter k = (*first_list).begin(); k != (*first_list).end(); ++k){
                    v1 = *k;
                    v2 = (++k == (*first_list).end()) ? (*first_list).front() : *k;
                    --k;
                    //get the neighbor of polygon i with edge v1--->v2
                    nbour = get_nbour(i, v1, v2);
                    //skip the edge if it has no neighbor
                    if(nbour == 0)
                        continue;
                    //calculate the color distance of this neighbor
                    color_dist = dist_cyl(i->c,nbour->c);
                    //update the nearest neighbor in color
                    if (color_dist < nearest_color_dist){
                        nearest_color_dist = color_dist;
                        nearest_nbour = nbour;
                    }
                }
                //merge polygon i with nearest polygon in color
                //remove both polygons and add the merged polygon
                small_polys.push_back(merge_utility(i, nearest_nbour,F));
                small_polys.remove(nearest_nbour);
            }
            small_polys.pop_front();
        }
    }

    // check if Two Edges are Connected or Not
    //Prepared By:Mahmoud Kareem
    bool Net::connected(polygon *P,Edge *e1, Edge*e2)
    {
        if ( e1->S != e2->S && e1->S != e2->E && e1->E != e2->S && e1->E != e2->E )
       {
            vertex *next=P->next_vertex(e1->S);
            vertex * prev=P->prev_vertex(e1->E);
            if ( e2->S !=next && e2->S !=prev &&e2->E !=next && e2->E !=prev  )
                  return false;  //not Connected
            else
                  return true;  // Connected
        }
       else
            return true;  //connected
    }


    // check if There is a child Polygon vertex  Exist in The Given Polygon P
    //Prepared By:Mahmoud Kareem
    bool Net::child_Polygon_Exist(std::list<std::list<vertex*>> vertices, polygon* P)
    {
        std::list<std::list<vertex*>>::iterator nested_iterator;  // Nested Iterator to iterate over Nested Lists
        nested_iterator= (vertices).begin(); //first List [ outer List ] not needed Here
        nested_iterator++;  //first Child if Exist
        for ( nested_iterator ; nested_iterator !=(vertices).end() ;nested_iterator ++)
        {
           std::list<vertex*>& list_ptr=*nested_iterator ; // Pointer to  List of The Vertices of The Polygon
           // loop over all vertices of The Polygon List
           for ( auto it =list_ptr.begin(); it!= list_ptr.end(); it++)
           {
               // for Each Vertex check if it is Exist in Polygon P or not
                  if ( P->PIP (  (*it)->V_point() )  )
                      {
                            // Vertex Exist in Polygon  return True and Exite
                            return true;
                      }
               }
        }
        return false; // There is no Child or no child Vertex Exist in The Given Polygon
    }

    //Prepared By:Mahmoud Kareem
    bool Net::edges_contracted(Edge *e1, Edge *e2,frame&F)
    {
         int x,y;
        // Calculate Sensitivity Region for e1
          diplacement D1=  measure_diplacement( e1->S , e1->E );
        // Vertex e1.S after diplacement by deltaX , deltaY
          x= int( e1->S->x +D1.deltaX );
          y= int(e1->S->y+D1.deltaY);
        // check if points is inside The Image Frame
        coordinates_check(&x,&y,F.height,F.width);
        vertex S1_disp( x,y ) ;

         // Vertex e1.E after diplacement by deltaX , deltaY
          x= int(e1->E->x +D1.deltaX);
          y= int(e1->E->y+D1.deltaY);
           // check if points is inside The Image Frame
          coordinates_check(&x,&y,F.height,F.width);
        vertex E1_disp( x,y ) ;

         vertex center1( (e1->S->x +e1->E->x)/2 ,  (e1->S->y +e1->E->y)/2  );
           x=int(center1.x +D1.deltaX);
          y=int(center1.y +D1.deltaY);
           // check if points is inside The Image Frame
         vertex center1_disp(x,y);
        // Calculate Sensitivity Region for e2
        diplacement D2=  measure_diplacement( e2->S , e2->E );
        // Vertex e2.S after diplacement by deltaX , deltaY
          x= int(e2->S->x +D2.deltaX);
          y= int(e2->S->y+D2.deltaY);
        //check if points is inside The Image Frame
        coordinates_check(&x,&y,F.height,F.width);
        vertex S2_disp( x,y ) ;

         // Vertex e2.E after diplacement by deltaX , deltaY
          x= int(e2->E->x +D2.deltaX);
          y= int(e2->E->y+D2.deltaY);
          coordinates_check(&x,&y,F.height,F.width);
        vertex E2_disp( x,y ) ;

        vertex center2( (e2->S->x +e2->E->x)/2 ,  (e2->S->y +e2->E->y)/2  );
          x=int(center2.x +D2.deltaX);
          y=int(center2.y +D2.deltaY);
           // check if points is inside The Image Frame
          vertex center2_disp(x,y);

        // e1 Sensetivity Region Polygon
        std::list<std::list<vertex*>> LLT;
        std::list<vertex*> LT;
        LT.push_back(e1->E );
        LT.push_back(e1->S );
        LT.push_back(&S1_disp );
        LT.push_back(&E1_disp);
        LLT.push_back( LT );
        polygon T1(LLT);

        // e2 Sensetivity Region Polygon
        LT.push_back(e2->E );
        LT.push_back(e2->S );
        LT.push_back(&S2_disp );
        LT.push_back(&E2_disp);
        LLT.push_back( LT );
        polygon T2(LLT);

        // check if The Two Sensitivity Regions are Overlapped or Not   [Condition may be Changed]
        //    if (  ( T.PIP(S2_disp.V_point() ) || T.PIP(E2_disp.V_point() ) )  && (TT.PIP(e1->S->V_point()) && T.PIP(e1->E->V_point()) )  )
           // if (  ( T.PIP(S2_disp.V_point() ) || T.PIP(E2_disp.V_point() ) )  || ( TT.PIP(S1_disp.V_point() )   || TT.PIP(E1_disp.V_point() )   )  )
        if (  ( T1.PIP(center2.V_point() ) || T1.PIP(center2_disp.V_point() ) )  && ( T2.PIP(center1.V_point() )   || T2.PIP(center1_disp.V_point() )   )  )
             return true;   //edges are Contracted
        else
            return false;  //edges not Contracted
    }

//    //split contracted polygons
//    //Prepared By:Mahmoud Kareem
//    void Net::split_polygons(frame F)
//    {
//        // loop over all Polygons in Polys List
//      std::list<polygon*>::iterator it;
//      static int m;
//      for ( it=polys.begin() ;it!= polys.end(); it++ )
//      {
//        // for Each Polygon Loop over all Child Polygons
//        polygon * Pi= *it;  // pointer to Each Polygon in The Poly list
//
//        std::list<std::list<vertex*> >::iterator nested_list_itr =  Pi->myVertices.begin() ; // First List [ outer Polygon List ]
//        nested_list_itr++;  // now we deal with first Child Polygons
//        // Loop over all  Child Polygons Lists of The Polygon Pi
//        static int q;
//        for (nested_list_itr ; nested_list_itr != Pi->myVertices.end();  nested_list_itr++)
//        {
//            //  [ Pointer to a Child Polygon List ]
//            std::list<vertex*> single_list_pointer = *nested_list_itr;
//            std::list<vertex*>::iterator itt; // Normal Iterator
//            // Loop over all Vertices of The Child Polygon
//               Edge e1,e2;
//               bool splitted=false;
//               polygon *PN, *PP1, *PP2;
//            for ( itt=single_list_pointer.begin(); itt != single_list_pointer.end(); itt++ )
//            {
//             std::list<vertex*>::iterator Temp;
//               // Construct  Edge e1
//               e1.S=*itt;  //one Vertex
//               Temp=itt;
//               Temp++;
//                if ( Temp == single_list_pointer.end() )
//                    e1.E= *single_list_pointer.begin();  //Next Vertex
//                else
//                   e1.E= *Temp;  //Next Vertex
//
//               // Search For Another Edge e2 make Construction with e1
//              std::list<vertex*>::iterator vv =itt;
//               for ( vv=single_list_pointer.begin(); vv != single_list_pointer.end(); vv++ )
//               {
//               // Construct  Edge e2
//               e2.S=*vv;  //one Vertex
//               Temp=vv;
//               Temp++;
//               if ( Temp == single_list_pointer.end() )
//                e2.E= *single_list_pointer.begin();  //Next Vertex
//                else
//                     e2.E=*Temp;
//                  //  Neigbour polygon of e1
//              //    if ( e1.E->x > e1.E->x)
//                // polygon *Pn1=get_nbour(Pi,e1.E,e1.S );
//                 //   Neigbour polygon of e2
//        //      polygon *Pn2=get_nbour(Pi,e2.E,e2.S );
//               polygon *Pn2,*Pn1;
//               if ( ( e1.E->x > e1.S->x ) || ( e1.E->x == e1.S->x && e1.E->y < e1.S->y) )
//                {
//                    Pn1=get_nbour(Pi,e1.S,e1.E );
//                }
//                else
//                {
//                    Pn1=get_nbour(Pi,e1.E,e1.S );
//                }
//                if ( ( e2.E->x > e2.S->x ) || ( e2.E->x == e2.S ->x&& e2.E->y < e2.S->y) )
//                {
//                    Pn2=get_nbour(Pi,e2.S,e2.E );
//                }
//                else
//                {
//                    Pn2=get_nbour(Pi,e2.E,e2.S );
//                }
//
//               if ( Pn1 == Pn2  &&  !connected( Pn2, &e1 , &e2 )  && edges_contracted ( &e1 , & e2,F)   )
//               {
//                   // now  e1 and e2 are contracted
//                   // Check That edges Satisfy This Condition
//                   // no any another child  vertex inside The Region between These two Edges
//                   // if This Conditions are Achieved Then The Polygon Must be Splitted into Two Smallers Polygons P1 and P2
//
//                 std::list<std::list<vertex*>> LLT;  //List of List of Temp Polygon represnt The Region between e1 and e2
//                 // This Polygon only have one List [ outer Polygon Vertices ]
//                  std::list<vertex*> LT;  // List of Outer Vertices [ Clockwise ordre ]
//                  LT.push_back(e1.E );
//                  LT.push_back(e1.S );
//                  LT.push_back(e2.E );
//                  LT.push_back(e2.S );
//                  LLT.push_back( LT );
//                // Construct Temp Polygon to deal with Child Polygons if Founded
//                  polygon T(LLT);
//                 //check if Edges Have Same Neghbour and no child Vertex are Exist between Edges
//                 static int m;
//                 if ( !child_Polygon_Exist(  Pn1->myVertices , &T ) )
//                 {
//                     std::list<std::list<vertex*>> LL1;  //List of List of Polygon P1
//                     std::list<std::list<vertex*>> LL2; //List of List of Polygon P2
//                     std::list<vertex*> L1;                      //For P1
//                     std::list<vertex*> L2;                      // for P2
//                      // now we must Split The Polygon [ Pn1 or Pn2 it is the same one ] into two smaller polygons P1 and P2
//                      //single_list_pointer is Point to The Splitted Polygon
//
//                     std::list<vertex*> ::iterator P_start;  // Iterator to The Start Vertex of The Polygon
//                     std::list<vertex*> ::iterator P_end;  // Iterator to The Start Vertex of The Polygon
//                     std::list<vertex*> ::iterator n_it;  // Normal Iterator
//                  // First Fill The Outer Frame of P1 [ Clockwise ]
//                    P_start= std::find (single_list_pointer.begin(), single_list_pointer.end(), e1.E);
//                    P_end = std::find (single_list_pointer.begin(), single_list_pointer.end(), e2.S);
//                     for ( n_it = P_start ; n_it != P_end ; n_it ++)
//                     {
//                         if ( n_it == single_list_pointer.end() ) // if The List Ended
//                            n_it = single_list_pointer.begin(); //Start Again
//                        //Put The Vertex in The L1 [ clockwise ordre ]
//                         L1.push_front( *n_it);
//                     }
//                     L1.push_front( *n_it);
//
//                    // Fill The Outer Frame of P2 [ Clockwise ]
//                    P_start= std::find (single_list_pointer.begin(), single_list_pointer.end(), e2.E);
//                    P_end = std::find (single_list_pointer.begin(), single_list_pointer.end(), e1.S);
//
//                     for ( n_it = P_start  ;  n_it != P_end ;  ++n_it)
//                     {
//                         if ( n_it == single_list_pointer.end() ) // if The List Ended
//                            n_it = single_list_pointer.begin(); //Start Again
//                        //Put The Vertex in The L1 [ clockwise ordre ]
//                        if ( n_it == P_end )
//                            break;
//                         L2.push_front( *n_it);
//
//                     }
//                     L2.push_front( *n_it);
//
//
//                     LL1.push_back( L1 );  // Add The Outer Frame List to The List of List of P1
//                     LL2.push_back( L2 );  // Add The Outer Frame List  to The List of List of P1
//
//                 // Construct The Polygons P1 and P2
//                   polygon* P1 = new polygon( LL1 );
//                   polygon* P2 = new polygon( LL2 );
//
//                     // check The Child Polygons of The Polygon [ Pn1 ] where it will added to P1 or P2
//                    std::list<std::list<vertex*>>::iterator child_iter =Pn1->myVertices .begin();
//                    child_iter++;  //refer to First Child of The Polygon
//                  //Loop over All Vertices List [ Child Polygons ] of The Polygon Pn1
//                  for (child_iter  ; child_iter !=Pn1->myVertices.end() ; child_iter++)
//                    {
//                         std::list<vertex*> list_ptr=*child_iter;  // Pointer to  List of The Vertices of The Polygon
//                       // check if The Child Polygon is to be placed in P1 or P2
//                        if (  P1->PIP (   (*(list_ptr.begin()))->V_point()  )  )
//                           {
//                              // Add This child Polygon to P1
//                            //add This List to P1_vertices List
//                            P1->myVertices.push_back(list_ptr);
//                          }
//                      else
//                         {
//                            // Add This child Polygon to P2
//                            //add This List to P2_vertices List
//                           P2->myVertices.push_back(list_ptr);
//                        }
//                   }
//
//                    PP1=P1;
//                    PP2=P2;
//                    PN=Pn1;
//   //indicate That The Polygon are Splitted and Stop Any Iterating or Looping Over it
//      splitted=true;
//                 }
//
//               } // edges are contracted and not connected
//             if (splitted)
//                  break;
//               }// search for e2
//
//             if (splitted)
//              {
//                  break;   //Stop Any Iterating or Looping Over The Splittited Polygon
//              }
//            } // Child polygon Construction Edges Search  Loop
//
//          //check if there is a splitted Polygon
//     if ( splitted )
//         {
//             // update The Affected Data Due to Splite
//            splitted=false;
//             //update Polys List
//            std::list<polygon*>::iterator iter = it;
//           iter++;
//            polys.insert(it,PP1);
//            polys.insert(it,PP2);
//            polys.remove(PN);
//                //for all vertices of The Splitted Polygon  PN Remove it from its Polygon List and Add PP1 or PP2
//                 std::list<std::list<vertex*> >::iterator nest_iter ;
//                 //Loop over all Child  Polygons
//                 for (nest_iter =( PN->myVertices ).begin() ; nest_iter != (PN->myVertices ).end() ; nest_iter ++)
//                 {
//                    std::list<vertex*>::iterator i ;
//                    for ( i = (*nest_iter).begin() ; i != (*nest_iter).end() ; i++ )
//                    {
//                     if ( PP1->PIP( (*i)->V_point()  ) )
//                         {
//                             // Add P1 before Pn1 [ which will be Deleted ]
//                           (*i)->add_polygon(PP1 , PN );
//                         (*i)->remove_polygon( PN );
//                         }
//                    else
//                        {
//                            // Add P2 before Pn1 [ which will be Deleted ]
//                         (*i)->add_polygon(PP2 , PN );
//                          (*i)->remove_polygon( PN);
//                        }
//                    }// Polygon Vertices Loop
//                 } //Polygons List Loop
//
//   //For The Original Polygon Pi remove The Splitted Polygon PN from it and add The Two Polygons PP1 and PP2 Lists [ Anti Clockwise ]
//    std::list<vertex*>LT1= *( PP1->myVertices.begin());
//    std::list<vertex*>LT2= *( PP2->myVertices.begin());
//    LT1.reverse();
//    LT2.reverse();
//    std::list<std::list<vertex*>>::iterator pp=Pi->myVertices.begin();
//       pp++;
//        Pi->myVertices.insert(pp,LT1 );
//        Pi->myVertices.insert( pp,LT2);
//      Pi->remove_poly_list(single_list_pointer);
//         }
//        }// Child Polygons Loop
//      } // Polygons Loop
//
//}


//  void Net::insert_polygons(frame F)
//    {
//        float num_of_pix = 0;
//        int pix_same_col = 0;
//        std::list<polygon*>::iterator it1 = polys.begin();
//        for (auto& i : polys)
//        {
//            auto j = i->myVertices.front();
//            //calculate square borders
//            int maxX = INT_MIN, maxY = INT_MIN, minX = INT_MAX, minY = INT_MAX;
//            for (auto k : j)
//            {
//                if (k->x > maxX)
//                    maxX = k->x;
//                if (k->x < minX)
//                    minX = k->x;
//                if (k->y > maxY)
//                    maxY = k->y;
//                if (k->y < minY)
//                    minY = k->y;
//            }
//            auto newlist = new std::list<vertex*>;
//            int num_ver = 0;
//            int iter = 10;
//            // calculate filling factor
//            for (int a = minY; a <= maxY ; a++)
//            {
//                for (int b = minX; b <= maxX ;b++)
//                {
//                    if (i->PIP(point(b, a)))
//                    {
//                        num_of_pix++;
//                        if (dist_cyl(color(F.data[b][a]), i->c) <= this->zeta_thresh)
//                        {
//                            pix_same_col++;
//                        }
//                        if (iter == 10)
//                        {
//                            num_ver++;
//                            auto newvertex = new vertex(b, a);
//                            newlist->emplace_back(newvertex);
//                            iter = 0;
//                        }
//                        iter++;
//                    }
//                }
//
//            }
//            if ((pix_same_col / num_of_pix) < fill_thresh)
//            {
//                //update polygon
//                auto L = new std::list<std::list<vertex*>>;
//                L->front() = *newlist;
//                polygon* pol = new polygon(*L);
//                pol->c = i->c;
//                pol->n = num_ver;
//                pol->area = i->area;
//                //copy children polys
//                int o = 0;
//                for (auto m : i->myVertices)
//                {
//                    if (o)
//                        pol->myVertices.emplace_back(m);
//                    o++;
//                }
//
//                it1 = polys.erase(it1);
//                polys.insert(it1, pol);
//                //update vertices
//
//                for (auto i : *newlist)
//                {
//                    verts.emplace_back(i);
//                }
//            }
//            else
//            {
//                delete newlist;
//            }
//
//            it1++;
//        }
//    }

    //delete unnecessary vertices
    //Prepared BY: Khaled M. Kamel
    void Net::delete_vertices(){
        typedef std::list<polygon*>::iterator poly_iter;
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;

        list_iter first_list;
        vertex *v1 , *v2 , *v3;
        polygon *e1_nbour, *e2_nbour;
        ver_iter temp;
        bool intersect;
        vertex *vs, *ve;
        float a1,b1,c1,a2,b2,c2;
        float s, t;
        float denom;

        //for all pairs of adjacent edges (e1) v1--->v2 , (e2) v2--->v3
        //remove v2 if : /**Aligned Edges**/
        //               e1 , e2 have angle difference < 5 degrees
        //               or
        //               /**Small Length Edge**/
        //               e1 has length < 10 pixels
        //               or
        //               /**Spike Edges**/
        //               e1 , e2 enclose an angle  < 10 degrees
        //do not remove v2 if : v2 is part of any other edge
        //                      or
        //                      v2 is part of a 3-vertex polygon (v1 v2 v3)

        //Aligned Edges threshold ( cos(5 degrees) )
        const float cos5 = 0.9961946981;
        //Small Length Edge threshold
        const int small_len = 10;
        //Spike Edges threshold
        const float cos10 = 0.984807753;

        double deltaA;
        point testpoint;
        //for all pairs of adjacent edges (e1) v1--->v2 , (e2) v2--->v3
        for(poly_iter i = polys.begin(); i != polys.end(); ++i){
            //iterate the outer edges only(thus covering all edges in the net)
            first_list = (*i)->myVertices.begin();
            //for polygons that have > 3 vertices
            if((*first_list).size() > 3){
                for(ver_iter j = (*first_list).begin(); j != (*first_list).end(); ++j){
                    //current vertex
                    v2 = *j;
                    //do not consider net corners
                    if( !( (v2->x==0&&v2->y==0) ||
                           (v2->x==image_width&&v2->y==0) ||
                           (v2->x==image_width&&v2->y==image_height) ||
                           (v2->x==0&&v2->y==image_height)) ){

                        temp = j;
                        //previous vertex
                        v1 = (--temp == (*first_list).end()) ? (*first_list).back() : *temp;
                        temp = j;
                        //next vertex
                        v3 = (++temp == (*first_list).end()) ? (*first_list).front() : *temp;
                        e1_nbour = get_nbour(*i, v1, v2);
                        e2_nbour = get_nbour(*i, v2, v3);
                        //if v2 is part of e1 and e2 only
                        //in polygons that both have > 3 vertices in the outer edges
                        if( e1_nbour == e2_nbour &&
                           (e1_nbour == 0|| e1_nbour->myVertices.front().size() > 3) ){

                            vector e1(v2->x-v1->x, v2->y-v1->y);
                            vector e2(v3->x-v2->x, v3->y-v2->y);
                            vector e1_reversed(v1->x-v2->x, v1->y-v2->y);
                            /**Aligned Edges**/
                            //if e1 , e2 have angle difference < 5 degrees
                            /**Small Length Edge**/
                            //or if e1 has length < 10 pixels
                            /**Spike Edges**/
                            //or if e1 , e2 enclose an angle < 10 degrees
                            if( (e1.dot(e2)/(e1.mag()*e2.mag())) > cos5 ||
                                e1.mag() < small_len ||
                                (e1_reversed.dot(e2)/(e1_reversed.mag()*e2.mag())) > cos10 ){

                                /*check if the edge resulting from deleting the current edge (v1-->v3) intersects any other
                                  edge in the polygon*/
                                intersect = false;
                                for(poly_iter k = polys.begin(); k != polys.end(); ++k){
                                    for(list_iter m = (*k)->myVertices.begin(); m != (*k)->myVertices.end(); ++m){
                                        for(ver_iter n = (*m).begin(); n!= (*m).end(); ++n){
                                            vs = *n;
                                            ve = (++n == (*m).end()) ? (*m).front() : *n;
                                            --n;
                                            if(vs != v1 && vs != v3 && ve != v1 && ve != v3){
                                                a1 = v1->x - v3->x;
                                                b1 = vs->x - ve->x;
                                                c1 = vs->x - v3->x;
                                                a2 = v1->y - v3->y;
                                                b2 = vs->y - ve->y;
                                                c2 = vs->y - v3->y;
                                                denom = a1*b2 - a2*b1;
                                                if(denom != 0){
                                                    s = (c1*b2 - c2*b1) / denom;
                                                    t = (a1*c2 - a2*c1) / denom;
                                                    if(s>=0 && s<=1 && t>=0 && t<=1 ){
                                                        intersect = true;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                        if(intersect)
                                            break;
                                    }
                                    if(intersect)
                                        break;
                                }
                                if(!intersect){
                                    //remove v2 from polygon i neighbor -if it exists-
                                    if(e1_nbour != 0)
                                        e1_nbour->remove_vertex(v2);
                                    //remove v2 from polygon i
                                    j = --(*first_list).erase(j);
                                    //update the areas and color averages of current and neighbor polygons
                                    //calculate the of the triangle v1->v2->v3
                                    deltaA = 0.5 * std::abs(v1->x*v2->y - v2->x*v1->y
                                                         +v2->x*v3->y - v3->x*v2->y
                                                         +v3->x*v1->y - v1->x*v3->y);
                                    if(deltaA != 0){
                                        testpoint.update_coords(v2->x, v2->y);
                                        if((*i)->PIP(testpoint)){
                                            (*i)->c.R = ((*i)->c.R * (*i)->area + e1_nbour->c.R * deltaA)/((*i)->area + deltaA);
                                            (*i)->c.G = ((*i)->c.G * (*i)->area + e1_nbour->c.G * deltaA)/((*i)->area + deltaA);
                                            (*i)->c.B = ((*i)->c.B * (*i)->area + e1_nbour->c.B * deltaA)/((*i)->area + deltaA);
                                            computeHSI((*i)->c);
                                            (*i)->area += deltaA;
                                            e1_nbour->area -= deltaA;

                                        }
                                        else{
                                            e1_nbour->c.R = ((*i)->c.R * deltaA + e1_nbour->c.R * e1_nbour->area)/(deltaA + e1_nbour->area);
                                            e1_nbour->c.G = ((*i)->c.G * deltaA + e1_nbour->c.G * e1_nbour->area)/(deltaA + e1_nbour->area);
                                            e1_nbour->c.B = ((*i)->c.B * deltaA + e1_nbour->c.B * e1_nbour->area)/(deltaA + e1_nbour->area);
                                            computeHSI(e1_nbour->c);
                                            (*i)->area -= deltaA;
                                            e1_nbour->area += deltaA;
                                        }
                                    }
                                    //finally remove v2 from the net
                                    remove_vertex(v2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //calculate the four disparities associated with an edge v1--->v2 in poly
    //Prepared BY: Khaled M.Kamel
    edge_disp Net::edge_disparity(polygon* poly, polygon* nbour, vertex* v1, vertex* v2, frame F){

        edge_disp result;
        point v1_point(v1->x, v1->y);
        point v2_point(v2->x, v2->y);
        //v1--->v2 mid point
        point mid( (v1->x+v2->x)/2, (v1->y+v2->y)/2 );

        //v1--->v2 vector
        vector e( v2->x-v1->x, v2->y-v1->y);
        //half v1--->v2 length
        int len = e.mag()/2;
        //edge unit vector
        e = e.unit_vec();
        //normal to edge unit vector
        vector v1v2_normal = e.perpendicular();
        //initialize current pixel to non-existent one
        point pixel(image_width, image_height);
        result.H_L1 = 0;
        //compute H_L1
        for(int j = 0; j <= S_width; ++j){
            for(int i = 0; i <= len; ++i){
                /*
                    moving a point by a vector might give the same point
                    skip a point if it is equal to the previous point
                */
                pixel = v1_point.move(v1v2_normal*j+ e*i);
                if(//nbour->PIP(pixel) && //if the point is in neighbor's region
                    pixel.x>=0 && pixel.x<image_width && pixel.y>=0 && pixel.y<image_height &&
                    //color distance between current pixel and neighbor above color theshold
                    (dist_cyl( F.data[pixel.x][pixel.y], nbour->c) > zeta_thresh) &&
                    //color distance between current pixel and poly below or equal color theshold
                    (dist_cyl( F.data[pixel.x][pixel.y], poly->c) <= zeta_thresh) ){

                    ++result.H_L1;
                }
            }
        }

        pixel.update_coords(image_width, image_height);
        result.H_R1 = 0;
        //compute H_R1
        for(int j = 0; j <= S_width; ++j){
            for(int i = 0; i <= len; ++i){
                pixel = mid.move(v1v2_normal*j+ e*i);
                if( //nbour->PIP(pixel) && //if the point is in neighbor's region
                    pixel.x>=0 && pixel.x<image_width && pixel.y>=0 && pixel.y<image_height &&
                    //color distance between current pixel and neighbor above color theshold
                    (dist_cyl( F.data[pixel.x][pixel.y], nbour->c) > zeta_thresh) &&
                    //color distance between current pixel and poly below or equal color theshold
                    (dist_cyl( F.data[pixel.x][pixel.y], poly->c) <= zeta_thresh) ){

                    ++result.H_R1;
                }
            }
        }

        //v2--->v1 vector
        vector e_reversed( v1->x-v2->x, v1->y-v2->y);
        //edge unit vector
        e_reversed = e_reversed.unit_vec();
        //normal to edge unit vector
        vector v2v1_normal = e_reversed.perpendicular();
        pixel.update_coords(image_width, image_height);
        result.H_L2 = 0;
        //compute H_L2
        for(int j = 0; j <= S_width; ++j){
            for(int i = 0; i <= len; ++i){
                pixel = v2_point.move(v2v1_normal*j+ e_reversed*i);
                if( //nbour->PIP(pixel) && //if the point is in neighbor's region
                    pixel.x>=0 && pixel.x<image_width && pixel.y>=0 && pixel.y<image_height &&
                    //color distance between current pixel and neighbor above color theshold
                    (dist_cyl( F.data[pixel.x][pixel.y], nbour->c) > zeta_thresh) &&
                    //color distance between current pixel and poly below or equal color theshold
                    (dist_cyl( F.data[pixel.x][pixel.y], poly->c) <= zeta_thresh) ){

                    ++result.H_L2;
                }
            }
        }


        pixel.update_coords(image_width, image_height);
        result.H_R2 = 0;
        //compute H_R2
        for(int j = 0; j <= S_width; ++j){
            for(int i = 0; i <= len; ++i){
                pixel = mid.move(v2v1_normal*j+ e_reversed*i);
                if( //nbour->PIP(pixel) && //if the point is in neighbor's region
                    pixel.x>=0 && pixel.x<image_width && pixel.y>=0 && pixel.y<image_height &&
                    //color distance between current pixel and neighbor above color theshold
                    (dist_cyl( F.data[pixel.x][pixel.y], nbour->c) > zeta_thresh) &&
                    //color distance between current pixel and poly below or equal color theshold
                    (dist_cyl( F.data[pixel.x][pixel.y], poly->c) <= zeta_thresh) ){

                    ++result.H_R2;
                }
            }
        }


        return result;
    }


    //break edges that have an overall force below the force threshold with a disparity measure above segma threshold
    //Prepared BY: Khaled M. Kamel
    void Net::insert_vertices(frame F){
        typedef std::list<polygon*>::iterator poly_iter;
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;

        list_iter first_list;
        vertex *v1, *v2;
        float len;
        edge_disp D;
        float v1_force, v2_force;
        polygon* nbour;
        //for all edges in all polygons
        for(poly_iter i = polys.begin(); i != polys.end(); ++i){
            first_list = (*i)->myVertices.begin();
            for(ver_iter j = (*first_list).begin(); j != (*first_list).end(); ++j){
                v1 = *j;
                v2 = (++j == (*first_list).end()) ? (*first_list).front() : *j;
                --j;
                //calculate the length of the edge v1--->v2
                len = sqrt(pow(v2->x - v1->x, 2) + pow(v2->y - v1->y, 2));
                nbour = get_nbour(*i, v1, v2);
                //skip the edge if it has no neighbor or if has a length less than 20
                if(nbour == 0 || len < 20)
                    continue;
                //get the disparities associated with edge v1--->v2
                D = edge_disparity(*i, nbour, v1, v2, F);
                //calculate the overall forces on v1 and v2
                v1_force = std::abs((D.H_L1-D.H_R2)/(2*len));
                v2_force = std::abs((D.H_R1-D.H_L2)/(2*len));
                //calculate the disparity threshold (semga * area of sensitivity region)
                float disp_thresh = segma_thresh * len/2 * S_width;
                //for an edge having an overall force below the force threshold with a disparity measure above segma threshold
                if(   ((v1_force<force_thresh)&&(D.H_L1>disp_thresh||D.H_R2>disp_thresh))
                    ||((v2_force<force_thresh)&&(D.H_R1>disp_thresh||D.H_L2>disp_thresh)) ){

                    std::list<polygon*> poly_list;
                    poly_list.push_front(nbour);
                    poly_list.push_front(*i);
                    //create a Midpoint vertex
                    vertex* v = new vertex((v1->x+v2->x)/2, (v1->y+v2->y)/2, poly_list);
                    //add the vertex to the Net
                    verts.push_back(v);
                    deformable.push_back(v);
                    //add the vertex to the neighboring polygon before v1
                    nbour->add_vertex(v, v1);
                    //add the vertex to this polygon before v2
                    j = (*first_list).insert(++j, v);
                }
            }
        }
    }

    //A full maintenance cycle
    //Prepared BY : Khaled M.Kamel
    void Net::maintenance(frame F){
        //calculate the area for all polygons in Net
        std::list<polygon*>::iterator i;
        for(i = polys.begin(); i != polys.end(); ++i){
            (*i)->area = (*i)->compute_area();
        }
        //Merge the neighboring polygons that have similar average colors
        color_merge(F);
        size_merge(F);
//        //Split the polygons that suffer from contraction
//        split_polygons(F);
//        //Insert new polygons if needed
//        insert_polygons(F);
        //Delete the unnecessary vertices
        delete_vertices();
        //Insert new vertices if needed,
        insert_vertices(F);
    }

    void Net::print_polys(){
        typedef std::list<polygon*>::iterator poly_iter;
        typedef std::list<std::list<vertex*>>::iterator list_iter;
        typedef std::list<vertex*>::iterator ver_iter;

        std::ofstream file;
        file.open("polys.txt");
        for(poly_iter i = polys.begin(); i != polys.end(); ++i){
            for(list_iter j = (*i)->myVertices.begin(); j != (*i)->myVertices.end(); ++j){
                file << "(" <<int((*j).front()->x) << "," << int((*j).front()->y) << ") ";
                for(ver_iter k = ++(*j).begin(); k != (*j).end(); ++k){
                    file << "(" << int((*k)->x) << "," << int((*k)->y) << ") ";
                    file << "(" << int((*k)->x) << "," << int((*k)->y) << ") ";
                }
                file << "(" <<int((*j).front()->x) << "," << int((*j).front()->y) << ")" << '\n';
            }
        }
        file.close();

    }

    void Net::print_verts(){
        typedef std::list<vertex*>::iterator ver_iter;

        std::ofstream file;
        file.open("verts.txt");
        for(ver_iter i = verts.begin(); i != verts.end(); ++i){
            file << "(" << int((*i)->x) << "," << int((*i)->y) << ") ";
        }
        file.close();
    }
