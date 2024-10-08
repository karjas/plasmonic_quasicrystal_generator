#ifndef SHAPE_H
#define SHAPE_H

#include <iostream> 
#include <vector>
#include <string>
#include <complex.h>
#undef I

#include <SFML/Graphics.hpp>


using namespace std;

class coord{

    public:
        coord(double x, double y);
        coord();
        ~coord(){}

        double x,y;

        double cross(coord b);
        double dot(coord b);
        double angle();

        coord operator + (coord const &c2){
            coord out;
            out.x = x + c2.x;
            out.y = y + c2.y;
            return out;
        };
        
        coord operator - (coord const &c2){
            coord out;
            out.x = x - c2.x;
            out.y = y - c2.y;
            return out;
        };
        
        coord operator * (double const f){
            coord out;
            out.x = x*f;
            out.y = y*f;
            return out;
        };
        
        coord operator / (double const f){
            coord out;
            out.x = x/f;
            out.y = y/f;
            return out;
        };

        coord unit();
        double len();

        coord c8();
        coord c4();
        coord cn(int N);

};

// E_fields
class E_fields{

    public:
        E_fields(int symm);
        ~E_fields(){};

        double calc_field_at(coord r, string irr);

        map<string,vector<complex<double>>*> get_irreps(){return irreps;};
        
        double get_k_mag();
        double get_k0x(){return kk[0][0];};
        double get_k0y(){return kk[0][1];};
        double get_mindist(){return mindist;};
        
        void set_mindist(double new_min){mindist = new_min;};

    private:
        map<string,vector<complex<double>>*> irreps;

        double kk[20][2];
        int Cn;
        double mindist = 0.5;

};

class shape{
    public:
        shape(vector<coord*> c, int orientation, double period);
        ~shape(){};
 
        virtual vector<shape*> subdivide(vector<coord*>* cc);
        vector<coord*> get_corners(){return corners;};

        virtual sf::ConvexShape gen_polygon();

        virtual coord optimize_point(E_fields* field, string irr, vector<shape*>* debug);

        int get_orientation(){return orientation;};

        bool get_mode(){return minimize;};
        
        void change_mode(){
            if(minimize){ 
                minimize = false;
            } else {
                minimize = true;
            }
        }



    protected:
        vector<coord*> corners;
        int orientation = 2;
        double period = 1;
        int scale = 5;
        bool minimize = true;

    private:
        string type = "S";

};


class Triangle: public shape{

    public:
        Triangle(vector<coord*> corners, int orientation, double p);
        ~Triangle(){}

        coord optimize_point(E_fields* field, string irr, vector<shape*>* debug);
        vector<shape*> subdivide(vector<coord*>* cc);
        sf::ConvexShape gen_polygon();
    private:
        string type = "T";

};

class Rhombic: public shape{

    public:
        Rhombic(vector<coord*> corners, int orientation, double p);
        ~Rhombic(){}

        coord optimize_point(E_fields* field, string irr, vector<shape*>* debug);
        vector<shape*> subdivide(vector<coord*>* cc);
        
        sf::ConvexShape gen_polygon();
    private:
        string type = "R";

};

class Segment: public shape{

    public:
        Segment(vector<coord*> corners, int orientation, double p);
        ~Segment(){}

        coord optimize_point(E_fields* field, string irr, vector<shape*>* debug);
        vector<shape*> subdivide(vector<coord*>* cc);
        
        sf::ConvexShape gen_polygon();

        double get_arc_max(){return arc_max;};
        void set_arc_max(double new_arc){arc_max = new_arc;};

    private:
        string type = "C";
        double arc_max = 1;

};

class Loop: public shape{

    public:
        Loop(vector<coord*> corners, int orientation, double p);
        ~Loop(){}

        vector<shape*> subdivide(vector<coord*>* cc);
        
        sf::ConvexShape gen_polygon();
        coord optimize_point(E_fields* field, string irr, vector<shape*>* debug);

        double get_arc_max(){return arc_max;};
        void set_arc_max(double new_arc){arc_max = new_arc;};

        bool is_inside_loop(coord c);

    private:
        string type = "C";
        double arc_max = 1.5;

};
#endif
