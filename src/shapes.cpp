#include <iostream>
#include <vector>
#include <math.h>
#include <string> 
#include <complex.h>
#undef I
#include <cfloat>

#include "shapes.hpp"

#include <SFML/Graphics.hpp>

using namespace std;

coord::coord(double x_, double y_){
    x = x_;
    y = y_;
}

coord::coord(){
    x = 0;
    y = 0;
}

double coord::len(){
    return sqrt(x*x + y*y);
}

double coord::angle(){
    return atan2(y,x);
}

coord coord::unit(){
    coord u;
    u.x = x;
    u.y = y;
    double l = u.len();
    return u/l;
}

double coord::cross(coord b){
    return x*b.y - y*b.x;
}

double coord::dot(coord b){
    return x*b.x + y*b.y;
}

coord coord::cn(int N){
    
    double th = 2*M_PI/N;
    coord u;
    u.x = cos(th)*x - sin(th)*y;
    u.y = cos(th)*y + sin(th)*x;
    return u;
}

coord coord::c8(){
    
    double th = M_PI/4;
    coord u;
    u.x = cos(th)*x - sin(th)*y;
    u.y = cos(th)*y + sin(th)*x;
    return u;
}
coord coord::c4(){
    
    double th = M_PI/2;
    coord u;
    u.x = cos(th)*x - sin(th)*y;
    u.y = cos(th)*y + sin(th)*x;
    return u;
}

// double[2] version of c8 rotation
void c8_rotation(double r0[2], double r[2]){
    double x,y;
    x = r0[0];
    y = r0[1];

    double th = M_PI/4;

    r[0] = cos(th)*x - sin(th)*y;
    r[1] = cos(th)*y + sin(th)*x;
}

void cn_rotation(double r0[2], double r[2], int N){
    double x,y;
    x = r0[0];
    y = r0[1];

    double th = 2*M_PI/N;

    r[0] = cos(th)*x + sin(th)*y;
    r[1] = cos(th)*y - sin(th)*x;
}

E_fields::E_fields(int symm){
    Cn = symm;

    double th = 2*M_PI/symm;
    /*
    cout << "th : " << th << endl;
    double px1, px2, py1, py2;
    px1 = 1;
    px2 = cos(th);
    py1 = 0;
    py2 = sin(th);
    
    double det = px1*py2 - px2*py1;

    cout << "det : " << det << endl;
    */
    double kl = 2*M_PI;
    (void)th;
    kk[0][1] = sin(0.0)*kl;
    kk[0][0] = cos(0.0)*kl;
    //kk[0][0] = py2/det*2*M_PI;
    //kk[0][1] = -py1/det*2*M_PI;

    cout << "Mag k0 : " << kk[0][1] << endl;

    for(int i = 0; i<symm;i++){
        cn_rotation(kk[i],kk[i+1],symm);
        cout << "k[" << i << "] : " << kk[i][0] << ", " << kk[i][1] << ", mag : " << sqrt(kk[i][0]*kk[i][0]+kk[i][1]*kk[i][1]  ) << endl;
    }
    if(symm == 6){
        complex<double> arrB[12] = {
            -sqrt(3)/3, sqrt(3)/6, sqrt(3)/6,
            -sqrt(3)/3, sqrt(3)/6, sqrt(3)/6,
            0, -0.5, 0.5,
            0, -0.5, 0.5
        };
        vector<complex<double>>* irrB = new vector<complex<double>>;
        irreps["B"] = irrB;

        for(int i = 0; i < 12; i++)irrB->push_back(arrB[i]);

        complex<double> arrA[12] = {
            0, sqrt(3)/2, sqrt(3)/2,
            0, -sqrt(3)/2, -sqrt(3)/2,
            1, 0.5, -0.5,
            -1, -0.5, 0.5
        };

        vector<complex<double>>* irrA = new vector<complex<double>>;
        irreps["A"] = irrA;
        for(int i = 0; i < 12; i++)irrA->push_back(arrA[i]);

        complex<double> arrE[12] = {
            0.25i,0.25i,0.25i,
            0.25i,0.25i,0.25i,
            0.25,0.25,0.25,
            0.25,0.25,0.25
        };

        vector<complex<double>>* irrE = new vector<complex<double>>;
        irreps["E"] = irrE;
        for(int i = 0; i < 12; i++)irrE->push_back(arrE[i]);

        complex<double> arrF[12] = {
            -sqrt(3)/12, -sqrt(3)/6 - 0.25i, sqrt(3)/6 - 0.25i,
            sqrt(3)/12, sqrt(3)/6 + 0.25i, -sqrt(3)/6 + 0.25i,
            sqrt(3)*0.25i, -0.25, -0.25,
            -sqrt(3)*0.25i, 0.25, 0.25
        };

        vector<complex<double>>* irrF = new vector<complex<double>>;
        irreps["F"] = irrF;
        for(int i = 0; i < 12; i++)irrF->push_back(arrF[i]);    
    }

    if(symm == 8){
        
        complex<double> arrB[16] = {
            -1, sqrt(2)/2, 0, -sqrt(2)/2,
            1, -sqrt(2)/2, 0, sqrt(2)/2,
            0, sqrt(2)/2, -1, sqrt(2)/2,
            0, -sqrt(2)/2, 1, -sqrt(2)/2
        };
        vector<complex<double>>* irrB = new vector<complex<double>>;
        irreps["B"] = irrB;
        for(int i = 0; i < 16; i++)irrB->push_back(arrB[i]);

        complex<double> arrA[16] = {
            1,0,-1,-sqrt(2),
            -1,0,1,sqrt(2),
            -1, -sqrt(2), -1, 0,
            1,sqrt(2),1,0};

        vector<complex<double>>* irrA = new vector<complex<double>>;
        irreps["A"] = irrA;
        for(int i = 0; i < 16; i++)irrA->push_back(arrA[i]);

        complex<double> arrE[16] = {
            -1i,0,1,1.0 - 1i,
            -1i,0,1,1.0 - 1i,
            1i,1.0 + 1i,1,0,
            1i,1.0 + 1i,1,0};

        vector<complex<double>>* irrE = new vector<complex<double>>;
        irreps["E"] = irrE;
        for(int i = 0; i < 16; i++)irrE->push_back(arrE[i]);
        complex<double> arrF[16] = {
            -1,0,-1,sqrt(2)*1i,
            1,0,1,-sqrt(2)*1i,
            1,-sqrt(2)*1i,-1,0,
            -1,sqrt(2)*1i,1,0};

        vector<complex<double>>* irrF = new vector<complex<double>>;
        irreps["F"] = irrF;
        for(int i = 0; i < 16; i++)irrF->push_back(arrF[i]);
        
        complex<double> arrG[16] = {
            1i,0,1,-1.0 - 1i,
            1i,0,1,-1.0 - 1i,
            -1i,-1.0 + 1i,1,0,
            -1i,-1.0 + 1i,1,0};
        
        vector<complex<double>>* irrG = new vector<complex<double>>;
        irreps["G"] = irrG;
        for(int i = 0; i < 16; i++)irrG->push_back(arrG[i]);
    }
    else if(symm == 10){
        complex<double> arrB[20] = {
            sqrt(50 - 10*sqrt(5))/10,  // 1
            sqrt(50 - 10*sqrt(5))/10,  // 2
            -(sqrt(5 - sqrt(5))*(5*sqrt(2) + 3*sqrt(10)))/20, // 3
            (sqrt(5 - sqrt(5))*(sqrt(10) + 5*sqrt(2)))/10,    //4
            -(sqrt(5 - sqrt(5))*(5*sqrt(2) + 3*sqrt(10)))/20, // 5
            sqrt(50 - 10*sqrt(5))/10,                           // 6
            sqrt(50 - 10*sqrt(5))/10,                           // 7
            -(sqrt(5 - sqrt(5))*(5*sqrt(2) + 3*sqrt(10)))/20,   // 8
            (sqrt(5 - sqrt(5))*(sqrt(10) + 5*sqrt(2)))/10,    //9
            -(sqrt(5 - sqrt(5))*(5*sqrt(2) + 3*sqrt(10)))/20,   // 10
            -sqrt(5)/2 - 0.5,
            0.5 + sqrt(5)/2,
            -1,
            0,
            1,
            -sqrt(5)/2 - 0.5,
            0.5 + sqrt(5)/2,
            -1,
            0,
            1
        };
        vector<complex<double>>* irrB = new vector<complex<double>>;
        irreps["B"] = irrB;
        for(int i = 0; i < 20; i++)irrB->push_back(arrB[i]);

        complex<double> arrA[20] = {
            sqrt(5 - sqrt(5))*(5*sqrt(2) + 3*sqrt(10))/20,
            sqrt(50 - 10*sqrt(5))/10,
            -sqrt(50 - 10*sqrt(5))/10,
            -sqrt(5 - sqrt(5))*(5*sqrt(2) + 3*sqrt(10))/20,
            -sqrt(5 - sqrt(5))*(sqrt(10) + 5*sqrt(2))/10,
            -sqrt(5 - sqrt(5))*(5*sqrt(2) + 3*sqrt(10))/20,
            -sqrt(50 - 10*sqrt(5))/10,
            sqrt(50 - 10*sqrt(5))/10,
            sqrt(5 - sqrt(5))*(5*sqrt(2) + 3*sqrt(10))/20,
            sqrt(5 - sqrt(5))*(sqrt(10) + 5*sqrt(2))/10,
            -1,
            -sqrt(5)/2 - 0.5,
            -sqrt(5)/2 - 0.5,
            -1,
            0,
            1,
            0.5 + sqrt(5)/2,
            0.5 + sqrt(5)/2,
            1,
            0
        };
        vector<complex<double>>* irrA = new vector<complex<double>>;
        irreps["A"] = irrA;
        for(int i = 0; i < 20; i++)irrA->push_back(arrA[i]);
    }
    if(symm == 12){
        
        complex<double> arrB[24] = {
            sqrt(3), -1, 0, 1,
            -sqrt(3), 2, -sqrt(3), 1,
            0, -1, sqrt(3), -2,
            -1, sqrt(3), -2, sqrt(3),
            -1, 0, 1, -sqrt(3),
            2, -sqrt(3), 1, 0,
        };
        vector<complex<double>>* irrB = new vector<complex<double>>;
        irreps["B"] = irrB;
        for(int i = 0; i < 24; i++)irrB->push_back(arrB[i]/(4*sqrt(3)));
    }
}



double E_fields::calc_field_at(coord r, string irrep){

    complex<double> Ex = 0;
    complex<double> Ey = 0;

    double pw = 0;

    for(int i = 0; i < Cn; i++){
        pw = (r.x*kk[i][0] + r.y*kk[i][1]); 
        Ex += irreps[irrep][0][i]* exp(1i*pw);
        Ey += irreps[irrep][0][i+Cn]* exp(1i*pw);
        //Ex += irreps[irrep][0][2*i]*exp(1i*pw);
        //Ey += irreps[irrep][0][2*i+1]*exp(1i*pw);
    }

    double out = real(Ex*conj(Ex)) + real(Ey*conj(Ey));

    return out;
}

double E_fields::get_k_mag(){
    return sqrt(kk[0][0]*kk[0][0] + kk[0][1]*kk[0][1]);
}

//shape::shape(vector<pair<double,double>> corners_, int orientation){
shape::shape(vector<coord*> corners_, int orient, double p){

    for(auto c : corners_)corners.push_back(c);
    orientation = orient;
    period = p;
    scale = 5;

}
sf::ConvexShape shape::gen_polygon(){
    sf::ConvexShape Polygon;
    return Polygon;
}


sf::ConvexShape Triangle::gen_polygon(){

    sf::ConvexShape Polygon;
    Polygon.setPointCount(corners.size());
    int i = 0;
    for(auto c : corners){
        Polygon.setPoint(i,sf::Vector2f(c->x*scale, c->y*scale));
        i++;
    }

    Polygon.setFillColor(sf::Color(251, 180, 175));
    Polygon.setOutlineThickness(1);
    Polygon.setOutlineColor(sf::Color(0, 0, 0));
    return Polygon;
}

sf::ConvexShape Rhombic::gen_polygon(){
    sf::ConvexShape Polygon;
    Polygon.setPointCount(corners.size());
    int i = 0;
    for(auto c : corners){
        Polygon.setPoint(i,sf::Vector2f(c->x*scale, c->y*scale));
        i++;
    }

    Polygon.setFillColor(sf::Color(179, 205, 227));
    Polygon.setOutlineThickness(1);
    Polygon.setOutlineColor(sf::Color(0, 0, 0));
    return Polygon;
}
sf::ConvexShape Segment::gen_polygon(){
    sf::ConvexShape Polygon;

    Polygon.setPointCount(corners.size());
    int i = 0;
    for(auto c : corners){
        Polygon.setPoint(i,sf::Vector2f(c->x*scale, c->y*scale));
        i++;
    }
    
    Polygon.setFillColor(sf::Color(179, 205, 227));
    Polygon.setOutlineThickness(1);
    Polygon.setOutlineColor(sf::Color(0, 0, 0));
    return Polygon;
}

//Triangle::Triangle(vector<pair<double,double>> corners, int orientation) : shape(corners, orientation){}
Triangle::Triangle(vector<coord*> corners, int orientation, double p) : shape(corners, orientation, p){}

Rhombic::Rhombic(vector<coord*> corners, int orientation, double p) : shape(corners, 1, p){(void)orientation;}

Segment::Segment(vector<coord*> corners, int orientation, double p) : shape(corners, 1, p){(void)orientation;}

Loop::Loop(vector<coord*> corners, int orientation, double p) : shape(corners, orientation, p){}

vector<shape*> shape::subdivide(vector<coord*>* corners){
    vector<shape*> sh;
    (void)corners;
    return sh;
}

pair<double,double> unit_vec(pair<double,double> p){
    double x = p.first;
    double y = p.second;
    
    double l = sqrt(x*x + y*y);
    pair<double,double> out;
    out.first = x/l;
    out.second= y/l;
    return out;
}


vector<shape*> Triangle::subdivide(vector<coord*>* cc){
    
    vector<shape*> out;

    coord *c1,*c2,*c3;
    c1 = corners[0];
    c2 = corners[1];
    c3 = corners[2];

    coord p1, p2, p3;
    p1 = *corners[0];
    p2 = *corners[1];
    p3 = *corners[2];

    double a0 = (p2 - p1).len();
    double a1 = a0/(1 + sqrt(2));
    
    // Points coming from subdivision
    coord p4 = p1 + (p2 - p1).unit()*a1*sqrt(2);
    coord* c4 = new coord(p4.x, p4.y);
    cc->push_back(c4);

    coord p5 = p2 + (p3 - p2).unit()*a1;
    coord* c5 = new coord(p5.x, p5.y);
    cc->push_back(c5);

    coord p6 = p3 + (p2 - p3).unit()*a1;
    coord* c6 = new coord(p6.x, p6.y);
    cc->push_back(c6);

    coord p7 = p1 + (p3 - p1).unit()*a1;
    coord* c7 = new coord(p7.x, p7.y);
    cc->push_back(c7);

    coord p8 = p1 + (p2*0.5 + p3*0.5 - p1).unit()*a1;
    coord* c8 = new coord(p8.x, p8.y);
    cc->push_back(c8);

    out.push_back(new Triangle({c8,c1,c4},orientation, period)); // Triangle A
    out.push_back(new Triangle({c8, c6, c5},-1*orientation, period)); // Triangle C
    out.push_back(new Triangle({c6, c3, c7}, orientation, period)); // Triangle D
    if(orientation == 1){
        out.push_back(new Rhombic({c8,c4, c2, c5}, 1, period)); // Rhombic B
        out.push_back(new Rhombic({c1, c8, c6, c7}, 1, period)); // Rhombic E
    } else {
        out.push_back(new Rhombic({c8,c5, c2, c4}, 1, period)); // Rhombic B
        out.push_back(new Rhombic({c1, c7, c6, c8}, 1, period)); // Rhombic E
    }

    
    return out;
}


vector<shape*> Rhombic::subdivide(vector<coord*>* cc){
     
    vector<shape*> out;
    coord *c1,*c2,*c3,*c4;
    c1 = corners[0];
    c2 = corners[1];
    c3 = corners[2];
    c4 = corners[3];

    coord p1, p2, p3, p4;
    p1 = *corners[0];
    p2 = *corners[1];
    p3 = *corners[2];
    p4 = *corners[3];

    double a0 = (p2 - p1).len();
    double a1 = a0/(1 + sqrt(2));

    coord p5 = p1 + (p2 - p1).unit()*a1;
    coord* c5 = new coord(p5.x, p5.y);
    cc->push_back(c5);

    coord p6 = p3 + (p2 - p3).unit()*a1;
    coord* c6 = new coord(p6.x, p6.y);
    cc->push_back(c6);

    coord p7 = p3 + (p4 - p3).unit()*a1;
    coord* c7 = new coord(p7.x, p7.y);
    cc->push_back(c7);
    
    coord p8 = p1 + (p4 - p1).unit()*a1;
    coord* c8 = new coord(p8.x, p8.y);
    cc->push_back(c8);

    coord p9 = p8 + (p2 - p1).unit()*a1;
    coord* c9 = new coord(p9.x, p9.y);
    cc->push_back(c9);
    
    coord p10 = p6 + (p1 - p2).unit()*a1;
    coord* c10 = new coord(p10.x, p10.y);
    cc->push_back(c10);

    out.push_back(new Rhombic({c1,c5,c9,c8},1, period)); // Rhombic A
    out.push_back(new Triangle({c9, c2, c5}, -1, period)); // Triangle B
    out.push_back(new Triangle({c10, c2, c6}, 1, period)); // Triangle C
    out.push_back(new Rhombic({c10,c6,c3,c7},1, period)); // Rhombic D
    out.push_back(new Triangle({c10, c4, c7}, -1, period)); // Triangle E
    out.push_back(new Triangle({c9, c4, c8}, 1, period)); // Triangle F
    out.push_back(new Rhombic({c2,c10,c4,c9},1, period)); // Rhombic G

    return out;
}

vector<shape*> Segment::subdivide(vector<coord*>* cc){
    vector<shape*> out;
    (void)cc;

    // Point of segment assumed at origin
    coord* c0 = corners[3];

    double th = c0->angle();
    double l = c0->len();
    int N = ceil(l/arc_max);

    double a = l/N;

    cout << N << " segments of size " << a << endl;

    double arc, dt;
    int M;

    for(int i = 0; i < N; i++){
        
        arc = th*(i+1)*a;
        M = ceil(arc/arc_max);
        dt = th/M;
        for(int j = 0; j < M; j++){

            coord* c1 = new coord(i*a*cos((j + 1)*dt), sin((j+1)*dt)*a*i);
            coord* c2 = new coord(i*a*cos(j*dt), i*a*sin(j*dt));
            coord* c3 = new coord((i+1)*a*cos(j*dt), sin(j*dt)*a*(i+1));
            coord* c4 = new coord((i+1)*a*cos((j+1)*dt), sin((j+1)*dt)*a*(i+1));

            cc->push_back(c1);
            cc->push_back(c2);
            cc->push_back(c3);
            cc->push_back(c4);

            out.push_back(new Segment({c1,c2,c3,c4},1, period)); 
        }
    }


    return out;
}

coord shape::optimize_point(E_fields* field, string irr, vector<shape*>* debug_shapes){
    /*
    coord u;
    (void)irr;
    (void)field;

    vector<coord*> scan_region;

    coord* c1 = corners[0];
    coord* c2 = corners[2];

    double min_dist = field->get_mindist();

    double l = c1->len();
    double dt = min_dist/l;

    double l2 = c2->len() - min_dist;
    double a1 = c1->angle() - dt;
    double a2 = c2->angle() + dt;

 
    scan_region.push_back(new coord(l*cos(a1), l*sin(a1)));
    scan_region.push_back(new coord(l*cos(a2), l*sin(a2)));
    scan_region.push_back(new coord(l2*cos(a2), l2*sin(a2)));
    scan_region.push_back(new coord(l2*cos(a1), l2*sin(a1)));

    debug_shapes->push_back(new Segment(scan_region, 1));

    coord uout;

    int N_steps = 25;
    double ddt = (a2 - a1)/N_steps;
    double ddl = (l2 - l)/N_steps;

    double t, ll;
    double th = 1;

    double ef = DBL_MAX;
    double ef2;

    for(int i = 0; i < N_steps - 1; i++){
        t = a1 + i*ddt; 
        for(int j = 0; j < N_steps - 1; j++){
            ll = l + j*ddl;
            u.x = cos(t)*ll;
            u.y = sin(t)*ll;
            ef2 = field->calc_field_at(u,irr);
            if(ef2 < ef){
                ef = ef2;
                uout = u;
            }
        }
    }

    if(ef > th){
        uout.x = DBL_MAX;
        uout.y = DBL_MAX;
    }

    return uout;

    */
    coord u;
    (void)irr;
    (void)field;

    int N = corners.size(); 

    vector<coord*> scan_region;
    for(auto c : corners){
        scan_region.push_back(new coord(c->x, c->y));
    }

    double shift = 0.1;
    (void)shift;

    coord c1, c2, cmid;
    coord norm, line, l1, l2;
    double nn, a1, a2;
    for(int i = 0; i < N; i++){
        c1 = *scan_region[i];
        
        c2 = *scan_region[(i+1)%N];


        cmid = c1*0.5 + c2*0.5;
        line = c2 - c1;
        norm = line.c4().unit()*orientation;
        nn = cmid.cross(norm);


        (void)c1;
        (void)c2;

        if(nn < 0){
            //l1 = (*corners[(i+N-1)%N] - *corners[i]).unit(); 
            l1 = (*scan_region[(i+N-1)%N] - *scan_region[i]).unit(); 
            //l2 = (*corners[(i+2)%N] - *corners[(i+1)%N]).unit();
            l2 = (*scan_region[(i+2)%N] - *scan_region[(i+1)%N]).unit();

            a1 = l1.dot(norm);
            a2 = l2.dot(norm);


            scan_region[i]->x += l1.x*shift/a1;
            scan_region[i]->y += l1.y*shift/a1;

            scan_region[(i+1)%N]->x += l2.x*shift/a2;
            scan_region[(i+1)%N]->y += l2.y*shift/a2;
        }
    }
    if(N == 3){
        debug_shapes->push_back(new Triangle(scan_region, orientation, period));
    }
    if(N == 4){
        debug_shapes->push_back(new Rhombic(scan_region, orientation, period));
    }

    int N_steps = 25;

    l1 = (*scan_region[1] - *scan_region[0]);
    double stepsize1 = l1.len()/N_steps;
    l1 = l1.unit()*stepsize1;

    l2 = (*scan_region[N-1] - *scan_region[0]);
    double stepsize2 = l2.len()/N_steps;
    l2 = l2.unit()*stepsize2;

    double ef = DBL_MAX;
    double ef2;

    coord uout;

    double th = 0.5;

    if(N == 3){
        for(int i = 0; i < N_steps - 1; i++){
            for(int j = 0; j < N_steps - i - 1; j++){
                u = l1*i + l2*j + *scan_region[0];
                ef2 = field->calc_field_at(u,irr);
                if(ef2 < ef){
                    ef = ef2;
                    uout = u;
                }
            }
        }
    }
    if(N == 4){
        for(int i = 0; i < N_steps - 1; i++){
            for(int j = 0; j < N_steps - 1; j++){
                u = l1*i + l2*j + *scan_region[0];
                ef2 = field->calc_field_at(u,irr);
                if(ef2 < ef){
                    ef = ef2;
                    uout = u;
                }
            }
        }
    }

    //TODO ADD further optimization..? Gradient descent etc

    if(ef > th) uout.x = DBL_MAX;

    return uout;
}


vector<shape*> Loop::subdivide(vector<coord*>* corners){
    vector<shape*> sh;
    cout << "Class (Loop) not subdiviseable" << endl;
    (void)corners;
    return sh;
}


sf::ConvexShape Loop::gen_polygon(){
    sf::ConvexShape Polygon;

    Polygon.setPointCount(corners.size());
    int i = 0;
    for(auto c : corners){
        Polygon.setPoint(i,sf::Vector2f(c->x*scale, c->y*scale));
        i++;
    }
    
    Polygon.setFillColor(sf::Color(179, 205, 0));
    Polygon.setOutlineThickness(1);
    Polygon.setOutlineColor(sf::Color(0, 0, 0));
    return Polygon;
}
coord Triangle::optimize_point(E_fields* field, string irr, vector<shape*>* debug_shapes){

    coord u;
    (void)irr;
    (void)field;

    int N = corners.size(); 

    vector<coord*> scan_region;
    for(auto c : corners){
        scan_region.push_back(new coord(c->x, c->y));
    }

    double shift = 0.1;
    (void)shift;

    coord c1, c2, cmid;
    coord norm, line, l1, l2;
    double nn, a1, a2;
    for(int i = 0; i < N; i++){
        c1 = *scan_region[i];
        
        c2 = *scan_region[(i+1)%N];


        cmid = c1*0.5 + c2*0.5;
        line = c2 - c1;
        norm = line.c4().unit()*orientation;
        nn = cmid.cross(norm);


        (void)c1;
        (void)c2;

        if(nn < 0){
            //l1 = (*corners[(i+N-1)%N] - *corners[i]).unit(); 
            l1 = (*scan_region[(i+N-1)%N] - *scan_region[i]).unit(); 
            //l2 = (*corners[(i+2)%N] - *corners[(i+1)%N]).unit();
            l2 = (*scan_region[(i+2)%N] - *scan_region[(i+1)%N]).unit();

            a1 = l1.dot(norm);
            a2 = l2.dot(norm);


            scan_region[i]->x += l1.x*shift/a1;
            scan_region[i]->y += l1.y*shift/a1;

            scan_region[(i+1)%N]->x += l2.x*shift/a2;
            scan_region[(i+1)%N]->y += l2.y*shift/a2;
        }
    }
    if(N == 3){
        debug_shapes->push_back(new Triangle(scan_region, orientation, period));
    }
    if(N == 4){
        debug_shapes->push_back(new Rhombic(scan_region, orientation, period));
    }

    int N_steps = 25;

    l1 = (*scan_region[1] - *scan_region[0]);
    double stepsize1 = l1.len()/N_steps;
    l1 = l1.unit()*stepsize1;

    l2 = (*scan_region[N-1] - *scan_region[0]);
    double stepsize2 = l2.len()/N_steps;
    l2 = l2.unit()*stepsize2;

    double ef = DBL_MAX;
    double ef2;

    coord uout;

    double th = 0.5;

    if(N == 3){
        for(int i = 0; i < N_steps - 1; i++){
            for(int j = 0; j < N_steps - i - 1; j++){
                u = l1*i + l2*j + *scan_region[0];
                ef2 = field->calc_field_at(u,irr);
                if(ef2 < ef){
                    ef = ef2;
                    uout = u;
                }
            }
        }
    }
    if(N == 4){
        for(int i = 0; i < N_steps - 1; i++){
            for(int j = 0; j < N_steps - 1; j++){
                u = l1*i + l2*j + *scan_region[0];
                ef2 = field->calc_field_at(u,irr);
                if(ef2 < ef){
                    ef = ef2;
                    uout = u;
                }
            }
        }
    }

    //TODO ADD further optimization..? Gradient descent etc

    if(ef > th) uout.x = DBL_MAX;

    return uout;
}
coord Rhombic::optimize_point(E_fields* field, string irr, vector<shape*>* debug_shapes){
    /*
    coord u;
    (void)irr;
    (void)field;

    vector<coord*> scan_region;

    coord* c1 = corners[0];
    coord* c2 = corners[2];

    double min_dist = field->get_mindist();

    double l = c1->len();
    double dt = min_dist/l;

    double l2 = c2->len() - min_dist;
    double a1 = c1->angle() - dt;
    double a2 = c2->angle() + dt;

 
    scan_region.push_back(new coord(l*cos(a1), l*sin(a1)));
    scan_region.push_back(new coord(l*cos(a2), l*sin(a2)));
    scan_region.push_back(new coord(l2*cos(a2), l2*sin(a2)));
    scan_region.push_back(new coord(l2*cos(a1), l2*sin(a1)));

    debug_shapes->push_back(new Segment(scan_region, 1));

    coord uout;

    int N_steps = 25;
    double ddt = (a2 - a1)/N_steps;
    double ddl = (l2 - l)/N_steps;

    double t, ll;
    double th = 1;

    double ef = DBL_MAX;
    double ef2;

    for(int i = 0; i < N_steps - 1; i++){
        t = a1 + i*ddt; 
        for(int j = 0; j < N_steps - 1; j++){
            ll = l + j*ddl;
            u.x = cos(t)*ll;
            u.y = sin(t)*ll;
            ef2 = field->calc_field_at(u,irr);
            if(ef2 < ef){
                ef = ef2;
                uout = u;
            }
        }
    }

    if(ef > th){
        uout.x = DBL_MAX;
        uout.y = DBL_MAX;
    }

    return uout;

    */
    coord u;
    (void)irr;
    (void)field;

    int N = corners.size(); 

    vector<coord*> scan_region;
    for(auto c : corners){
        scan_region.push_back(new coord(c->x, c->y));
    }

    double shift = 0.1;
    (void)shift;

    coord c1, c2, cmid;
    coord norm, line, l1, l2;
    double nn, a1, a2;
    for(int i = 0; i < N; i++){
        c1 = *scan_region[i];
        
        c2 = *scan_region[(i+1)%N];


        cmid = c1*0.5 + c2*0.5;
        line = c2 - c1;
        norm = line.c4().unit()*orientation;
        nn = cmid.cross(norm);


        (void)c1;
        (void)c2;

        if(nn < 0){
            //l1 = (*corners[(i+N-1)%N] - *corners[i]).unit(); 
            l1 = (*scan_region[(i+N-1)%N] - *scan_region[i]).unit(); 
            //l2 = (*corners[(i+2)%N] - *corners[(i+1)%N]).unit();
            l2 = (*scan_region[(i+2)%N] - *scan_region[(i+1)%N]).unit();

            a1 = l1.dot(norm);
            a2 = l2.dot(norm);


            scan_region[i]->x += l1.x*shift/a1;
            scan_region[i]->y += l1.y*shift/a1;

            scan_region[(i+1)%N]->x += l2.x*shift/a2;
            scan_region[(i+1)%N]->y += l2.y*shift/a2;
        }
    }
    if(N == 3){
        debug_shapes->push_back(new Triangle(scan_region, orientation, period));
    }
    if(N == 4){
        debug_shapes->push_back(new Rhombic(scan_region, orientation, period));
    }

    int N_steps = 25;

    l1 = (*scan_region[1] - *scan_region[0]);
    double stepsize1 = l1.len()/N_steps;
    l1 = l1.unit()*stepsize1;

    l2 = (*scan_region[N-1] - *scan_region[0]);
    double stepsize2 = l2.len()/N_steps;
    l2 = l2.unit()*stepsize2;

    double ef = DBL_MAX;
    double ef2;

    coord uout;

    double th = 0.5;

    if(N == 3){
        for(int i = 0; i < N_steps - 1; i++){
            for(int j = 0; j < N_steps - i - 1; j++){
                u = l1*i + l2*j + *scan_region[0];
                ef2 = field->calc_field_at(u,irr);
                if(ef2 < ef){
                    ef = ef2;
                    uout = u;
                }
            }
        }
    }
    if(N == 4){
        for(int i = 0; i < N_steps - 1; i++){
            for(int j = 0; j < N_steps - 1; j++){
                u = l1*i + l2*j + *scan_region[0];
                ef2 = field->calc_field_at(u,irr);
                if(ef2 < ef){
                    ef = ef2;
                    uout = u;
                }
            }
        }
    }

    //TODO ADD further optimization..? Gradient descent etc

    if(ef > th) uout.x = DBL_MAX;

    return uout;
}
coord Segment::optimize_point(E_fields* field, string irr, vector<shape*>* debug_shapes){

    coord u;
    (void)irr;
    (void)field;

    int N = corners.size(); 

    vector<coord*> scan_region;
    for(auto c : corners){
        scan_region.push_back(new coord(c->x, c->y));
    }

    double shift = 0.1;
    (void)shift;

    coord c1, c2, cmid;
    coord norm, line, l1, l2;
    double nn, a1, a2;
    for(int i = 0; i < N; i++){
        c1 = *scan_region[i];
        
        c2 = *scan_region[(i+1)%N];


        cmid = c1*0.5 + c2*0.5;
        line = c2 - c1;
        norm = line.c4().unit()*orientation;
        nn = cmid.cross(norm);


        (void)c1;
        (void)c2;

        if(nn < 0){
            //l1 = (*corners[(i+N-1)%N] - *corners[i]).unit(); 
            l1 = (*scan_region[(i+N-1)%N] - *scan_region[i]).unit(); 
            //l2 = (*corners[(i+2)%N] - *corners[(i+1)%N]).unit();
            l2 = (*scan_region[(i+2)%N] - *scan_region[(i+1)%N]).unit();

            a1 = l1.dot(norm);
            a2 = l2.dot(norm);


            scan_region[i]->x += l1.x*shift/a1;
            scan_region[i]->y += l1.y*shift/a1;

            scan_region[(i+1)%N]->x += l2.x*shift/a2;
            scan_region[(i+1)%N]->y += l2.y*shift/a2;
        }
    }
    if(N == 3){
        debug_shapes->push_back(new Triangle(scan_region, orientation, period));
    }
    if(N == 4){
        debug_shapes->push_back(new Rhombic(scan_region, orientation, period));
    }

    int N_steps = 25;

    l1 = (*scan_region[1] - *scan_region[0]);
    double stepsize1 = l1.len()/N_steps;
    l1 = l1.unit()*stepsize1;

    l2 = (*scan_region[N-1] - *scan_region[0]);
    double stepsize2 = l2.len()/N_steps;
    l2 = l2.unit()*stepsize2;

    double ef;
    if(minimize){
        ef = DBL_MAX;
    } else {
        ef = -DBL_MAX;
    }
    double ef2;

    coord uout;

    double th = 0.3;

    if(N == 3){
        for(int i = 0; i < N_steps - 1; i++){
            for(int j = 0; j < N_steps - i - 1; j++){
                u = l1*i + l2*j + *scan_region[0];
                ef2 = field->calc_field_at(u,irr);
                if(minimize){
                    if(ef2 < ef){
                        ef = ef2;
                        uout = u;
                    }
                } else {
                    if(ef2 > ef){
                        ef = ef2;
                        uout = u;
                    }
                }
            }
        }
    }
    if(N == 4){
        for(int i = 0; i < N_steps - 1; i++){
            for(int j = 0; j < N_steps - 1; j++){
                u = l1*i + l2*j + *scan_region[0];
                ef2 = field->calc_field_at(u,irr);
                if(minimize){
                    if(ef2 < ef){
                        ef = ef2;
                        uout = u;
                    }
                } else {
                    if(ef2 > ef){
                        ef = ef2;
                        uout = u;
                    }
                }
            }
        }
    }

    //TODO ADD further optimization..? Gradient descent etc
    

    if(minimize && (ef > th)) uout.x = DBL_MAX;
    if(!minimize && (ef < 30)) uout.x = DBL_MAX;

    return uout;
}

coord Loop::optimize_point(E_fields* field, string irr, vector<shape*>* debug_shapes){
    coord u;
    (void)irr;
    (void)field;

    int N = corners.size(); 
    (void)N;

    int N_steps = 25;
    (void)N_steps;
    (void)debug_shapes;

    double minx = DBL_MAX;
    double miny = DBL_MAX;
    double maxx = -DBL_MAX;
    double maxy = -DBL_MAX;

    for(auto c : corners){
        minx = min(c->x, minx);
        miny = min(c->y, miny);
        maxx = max(c->x, maxx);
        maxy = max(c->y, maxy);
    }
    double dx = (maxx - minx)/(N_steps - 1);
    double dy = (maxy - miny)/(N_steps - 1);
    coord uout;
    double th = DBL_MAX;
    double ef = DBL_MAX;
    double ef2;

    for(int i = 0; i < N_steps - 1; i++){
        for(int j = 0; j < N_steps - 1; j++){
            u.x = minx + i*dx;
            u.y = miny + i*dy;
            if(!is_inside_loop(u)) continue;

            ef2 = field->calc_field_at(u,irr);
            if(ef2 < ef){
                ef = ef2;
                uout = u;
            }
        }
    }
    if(ef > th) uout.x = DBL_MAX;

    return uout;
}


bool Loop::is_inside_loop(coord c){
 
    int N = corners.size();
    coord c1;
    coord c2;
    for(int i = 0; i < N; i++){
        c1 = *corners[i];
        c2 = *corners[(i+1)%N];
        //cout << c.x << ", " << c.y << "; "<<c1.x << ", " << c1.y << "; " << c2.x << ", " << c2.y << ", " << endl; 
        if((c.x - c1.x)*(c2.y - c1.y) + orientation*(c.y - c1.y)*(c2.x - c1.x) < 0) return false;       
    }
    return true;
}



