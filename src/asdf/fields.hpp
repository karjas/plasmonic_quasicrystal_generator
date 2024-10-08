#ifndef FIELDS_H
#define FIELDS_H

#include <vector>
#include <complex.h>
#undef I
#include <string>
#include <map>

#include "shapes.hpp"

class coord;

class E_fields{

    public:
        E_fields();
        ~E_fields(){};

        double calc_field_at(coord r, string irr);

 //       double calc_field_at_points(string fnam, string irrep);

//        double calc_coverage_r(string fnam, string irrep, double r, int N_points = 10);

    private:
        map<string,vector<complex<double>>*> irreps;

        double kk[8][2];

};

#endif
