#include <iostream>
#include "window.hpp"
#include "shapes.hpp"
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>


using namespace std;

int load_loops(string fnam, vector<shape*>* vsh, double p){
    
    (void)fnam;
    (void)vsh;

    ifstream file(fnam);
    string line;
    
    int loop_idx, curr_idx;
    double x,y;
    int scale = 1;
    
    int dx = 0;
    int dy = 0;
    
    if(file.is_open()){
        curr_idx = 0;
        (void) curr_idx;
        //vector<coord*>* cc = new vector<coord*>;
        vector<vector<coord*>*> cc;
        cc.push_back(new vector<coord*>);
        while(getline(file,line)){
            stringstream ss(line);
            ss >> loop_idx;
            ss.ignore();
            ss >> x;
            ss.ignore();
            ss >> y;

            if(curr_idx == loop_idx){
                cc[curr_idx]->push_back(new coord((x + dx)*scale,(y+dy)*scale));
            } else {
                Loop* sh = new Loop(*(cc[curr_idx]), -1, p);
                vsh->push_back(sh);
                cc.push_back(new vector<coord*>);
                cc[loop_idx]->push_back(new coord((x + dx)*scale,(y + dy)*scale));
                curr_idx += 1;
            }
        }   
    }  
    cout << "VSH size : " << vsh->size() << endl;

    return 0;
}

int main(){

    cout << "Hello World!" << endl;
   
    double x0 = 0;
    double y0 = 0;
    double b = 0;
    double p = 580e-9;
    //double a = pow(1+sqrt(2),6)*2.2;
    double a = 250;

    int symm = 12;
    double dt = M_PI/symm;
    dt = 0;
    double th = M_PI/symm;
    cout << a << endl;
    coord* c1 = new coord(x0 + cos(dt + th)*b,y0 + sin(dt + th)*b); 
    coord* c2 = new coord(x0 + cos(dt)*b,y0 + sin(dt)*b); 
    coord* c3 = new coord(x0 + cos(dt)*a,y0 + sin(dt)*a); 
    coord* c4 = new coord(x0 + cos(dt + th)*a,y0+ sin(dt + th)*a); 
    
    //coord* c1 = new coord(x0,y0); 
    //coord* c2 = new coord(x0 + a,y0); 
    //coord* c3 = new coord(x0 + a + cos(M_PI/4)*a,y0+ sin(M_PI/4)*a); 
    //coord* c4 = new coord(x0 + cos(M_PI/4)*a,y0+ sin(M_PI/4)*a); 


    //Rhombic* sh = new Rhombic(cc, 1, p);

    vector<coord*> cc = {c1,c2,c3,c4};

    Segment* sh = new Segment(cc, 1, p);
    //Loop* sh = new Loop(cc, 1, p);

    for(auto c : sh->get_corners()){
        cout << c->x << ", " << c->y << endl;
    }
    //vector<coord*>*cc2 = new vector<coord*>;

    vector<shape*> vsh = {sh};
    //vector<shape*> vsh = sh->subdivide(cc2);

    //vector<shape*> vsh2= {vsh[0],vsh[1],vsh[5],vsh[6]};

    //(void)vsh2;

    //mainWindow w = mainWindow(vsh2, symm, p);
    /*

    vector<shape*> vsh;

    (void)vsh;
    double p = 580e-9;
    int symm = 8;
    load_loops("loops.csv", &vsh, p);
    cout << vsh.size() << " loops loaded"  << endl;
    */
    mainWindow w = mainWindow(vsh, symm, p);

    return 0;
}
