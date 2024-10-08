#include "fields.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <math.h>
#include <complex.h>
#undef I


void c8_rotation(double r0[2], double r[2]){
    double x,y;
    x = r0[0];
    y = r0[1];

    double th = M_PI/4;

    r[0] = cos(th)*x - sin(th)*y;
    r[1] = cos(th)*y + sin(th)*x;
}


E_fields::E_fields(){

    double th = M_PI/4;

    double px1, px2, py1, py2;
    px1 = 1;
    px2 = cos(th);
    py1 = 0;
    py2 = sin(th);

    double det = px1*py2 - px2*py1;

    kk[0][0] = py2/det*2*M_PI;
    kk[0][1] = -px2/det*2*M_PI;

    for(int i = 0; i<7;i++){
        c8_rotation(kk[i],kk[i+1]);
    }


    complex<double> arrB[16] = {
        1,0,-1,sqrt(2),
        -1,0,1,-sqrt(2),
        -1,sqrt(2),-1,0,
        1,-sqrt(2),1,0};
    
    vector<complex<double>>* irrB = new vector<complex<double>>;
    irreps["B"] = irrB;
    for(int i = 0; i < 16; i++)irrB->push_back(arrB[i]);

    complex<double> arrA[16] = {
        1,0,-1,-sqrt(2),
        -1,0,1,sqrt(2),
        -1,-1,-1,0,
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
        1i,0,1,-1.0 - 1i,
        1i,0,1,-1.0 - 1i,
        -1i,-1.0 + 1i,1,0,
        -1i,-1.0 + 1i,1,0};

    vector<complex<double>>* irrF = new vector<complex<double>>;
    irreps["F"] = irrF;
    for(int i = 0; i < 16; i++)irrF->push_back(arrF[i]);
    
    complex<double> arrG[16] = {
        -1,0,-1,sqrt(2)*1i,
        1,0,1,-sqrt(2)*1i,
        1,-sqrt(2)*1i,-1,0,
        -1,sqrt(2)*1i,1,0};

    vector<complex<double>>* irrG = new vector<complex<double>>;
    irreps["G"] = irrG;
    for(int i = 0; i < 16; i++)irrG->push_back(arrG[i]);
}



double E_fields::calc_field_at(coord r, string irrep){

    complex<double> Ex = 0;
    complex<double> Ey = 0;

    double pw = 0;

    for(int i = 0; i < 8; i++){
        pw = (r.x*kk[i][0] + r.y*kk[i][1]); 
        Ex += irreps[irrep][0][i]* exp(1i*pw);
        Ey += irreps[irrep][0][i+8]* exp(1i*pw);
    }

    double out = real(Ex*conj(Ex)) + real(Ey*conj(Ey));

    return out;
}
