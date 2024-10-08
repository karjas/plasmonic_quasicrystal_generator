#include "window.hpp"
#include <string>
#include <fstream>
#include <cfloat>

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>


using namespace std;

mainWindow::mainWindow(vector<shape*> sh, int symm, double p){

    Cn = symm;
    period = p;
    
    w_ = new sf::RenderWindow(VideoMode(1200,800), "Field optimizer");

    points = new vector<coord*>;
    optimized_points = new vector<coord*>;
    debug_shapes = new vector<shape*>;

    fields = new E_fields(symm);

    bool found = false;
    for(auto s : sh){
        shapes.push_back(s);
        for(auto c : s->get_corners()){
            found = false;
            for(auto c2 : *points){
                if(c == c2) found = true;
            }
            if(!found) points->push_back(c);
        }
    }
    cout << "Initializing window with " << shapes.size() << " shapes!" << endl;
    draw();
    run();
}

void mainWindow::run(){
    while(w_->isOpen()){

        Event e;

        while(w_->pollEvent(e)){
            if(e.type == sf::Event::KeyPressed){
                if(sf::Keyboard::isKeyPressed(sf::Keyboard::S)){
                    cout << "Subdividing" << endl;
                    subdivide();
                    cout << "Now " << n_points() << " points" << endl;
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::W)){
                    cout << "Saving" << endl;
                    save_points();
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::Q)){
                    w_->close();
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::A)){
                    optimize_points("A");
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::B)){
                    optimize_points("B");
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::E)){
                    optimize_points("E");
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::F)){
                    optimize_points("F");
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::G)){
                    optimize_points("G");
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::M)){
                    change_mode();
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::D)){
                    if(drawing_points){ drawing_points = false;}
                    else{drawing_points = true;};
                }
                else if(sf::Keyboard::isKeyPressed(sf::Keyboard::T)){
                    if(drawing_points){ drawing_debug = false;}
                    else{drawing_debug = true;};
                }
            }
        }
        draw();
    }
    return;
}

void mainWindow::draw(){

    w_->clear(sf::Color::Green);
    for(auto s:shapes){
        w_->draw(s->gen_polygon()); 
    }

    if(drawing_debug){
        cout << "Drawing debug : " << debug_shapes->size() << endl;
        for(auto s:*debug_shapes){
            w_->draw(s->gen_polygon()); 
        }

    }
    if(drawing_points) draw_points();

    w_->display();
}

void mainWindow::subdivide(){

    vector<shape*> vsh;

    for(auto s : shapes){
        for(auto s2:s->subdivide(points)){
            vsh.push_back(s2);
        }
    }

    shapes.clear();
    for(auto s2 : vsh) shapes.push_back(s2);
    shape* sh0 = shapes[0];
    coord c1 = *(sh0->get_corners()[1]);
    coord c2 = *(sh0->get_corners()[2]);

    cout << "Distance : " << (c1-c2).len() << endl;

    return;
}


void mainWindow::print_points(){
    
    for(auto c : *points){
        cout << c->x << ", " << c->y << endl;
    }
    return;
}

double mainWindow::get_arc_max(){
    Segment *s = (Segment*)shapes[0];
    return s->get_arc_max();
}

bool mainWindow::get_mode(){
    bool mode = (Segment*)shapes[0]->get_mode();
    return mode;
}

void mainWindow::change_mode(){
    string mode1;
    string mode2;
    if(get_mode()){
        mode1 = "minimize";
        mode2 = "maximize";
    } else {
        mode1 = "maximize";
        mode2 = "minimize";
    }
    cout << "Changing mode from " << mode1 << " to " << mode2 << endl;
    for(auto s : shapes) s->change_mode(); 
    return;
}

void mainWindow::save_points(){

    string fnam;//; = last_opt + "_" + "points.csv";

    if(get_mode()){
        fnam = last_opt + "_min_points.csv";
    } else {
        fnam = last_opt + "_max_points.csv";
    }

    ofstream f(fnam);

    f << "Symm : " << get_symm() << endl;
    f << "kx : " << get_E_field()->get_k0x() << endl;
    f << "ky : " << get_E_field()->get_k0y() << endl;
    f << "kmag : " << get_E_field()->get_k_mag() << endl;
    f << "mindist  : " << get_E_field()->get_mindist() << endl;
    f << "Period : " <<  get_period()*1e9 << endl;

    f << "arc_max : " << get_arc_max() << endl;

    f << "PARTICLE POSITIONS" << endl;
    coord p1;
    cout << last_opt << endl;
    if(last_opt != "E"){
        f << "0 , 0" << endl;
    }
    //Cn = 1; // ONLY FOR LOADED SYSTEMS:
    bool on_axis = false;
    bool on_diag = false;
    double axis_dist = 0.2;
    double r;
    double dt = 2.0/Cn;
    double th0, th;
    if(last_opt == "N"){
        for(auto p : *points){
            if(p->y != 0){
                p1 = *p;
                if(p->x < 300){
                    for(int i = 0; i < Cn; i++){
                        f << p1.x*period<< ", " << p1.y*period << endl;
                        p1 = p1.cn(Cn);
                    }  
                }
            }
        }
    } else {
        for(auto p : *optimized_points){
            if(p->y != 0){
                p1 = *p;
                r = sqrt(p1.y*p1.y + p1.x*p1.x);        
                if(abs(p1.y) < axis_dist){
                    p1.y = 0;
                    on_axis = true;
                    on_diag = false;
                } else if(sin(M_PI/8)*p1.x - cos(M_PI/8)*p1.y - axis_dist < 0) {
                    p1.x = r*cos(M_PI/8);
                    p1.y = r*sin(M_PI/8);
                    on_axis = false;
                    on_diag = true;
                } else {
                    on_axis = false;
                    on_diag = false;
                }

                if(r > 300) continue;
                th0 = atan2(p1.x, p1.y)/M_PI;
                for(int i = 0; i < Cn; i++){
                    th = th0 + i*dt;
                    if(th > 1) th -= 2;
                    f << r*period << ", " << th << endl;
                    if(!on_axis && !on_diag){
                        th = -th0 + i*dt;
                        if(th > 1) th -= 2;
                        f << r*period << ", " << -th0 + i*dt << endl;
                    }
                }
            }
        }
    }

    f.close();

    return;
}



void mainWindow::optimize_points(string irr){
    
    cout << "Starting optimization" << endl;
    coord u;
    last_opt = irr;
    int N = shapes.size();
    int w = 20;
    int i = 0;
    for(auto s : shapes){
        coord* uu = new coord;
        u = s->optimize_point(fields, irr, debug_shapes);
        i++;
        if(u.x < 300){
            uu->x = u.x;
            uu->y = u.y;
            optimized_points->push_back(uu);

            if(i%(N/w) == 0){
                cout << i*100.0/N << "%" << endl;
            }
        }
    }
    cout << "Optimization complete!" << endl;

    return;
}

void mainWindow::draw_points(){

    cout << "Drawing points" << endl;
    cout << "Number of poitns : " << optimized_points->size() << endl;
    for(auto o: *optimized_points){
        sf::CircleShape shape(4);
        shape.setFillColor(sf::Color(150, 50, 250));
        shape.setOrigin(2,2);
        shape.setPosition(o->x*5, o->y*5);
        w_->draw(shape);
    } 

    w_->display();

    return;
}


