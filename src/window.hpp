#ifndef WINDOW_H
#define WINDOW_H

#include <iostream>
#include <string>
#include <vector>

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include "shapes.hpp"

using namespace std;
using namespace sf;

class mainWindow{

    public:
        mainWindow(vector<shape*>, int symm, double p);
        ~mainWindow(){}
        
        RenderWindow* getWindow(){return w_;};
        void add_shape(shape* s){shapes.push_back(s);};

        void subdivide();
        vector<coord*>* get_points(){return points;};
        int n_points(){return points->size();};
        void print_points();

        void save_points();
        E_fields* get_E_field(){return fields;};

        void optimize_points(string irr);
        vector<coord*>* get_optimized_points(){return optimized_points;};

        int get_symm(){return Cn;};

        void draw_points();

        double get_arc_max();
        double get_period(){return period;};

        bool get_mode();
        void change_mode();

    private:
        RenderWindow* w_;
    
        vector<shape*> shapes;
        vector<shape*>* debug_shapes;
        void draw();
        void run();
        vector<coord*>* points;
        vector<coord*>* optimized_points;
        E_fields* fields;
        string last_opt = "N";
        bool drawing_points = false;
        bool drawing_debug = false;
        int Cn;
        double period;
};




#endif
