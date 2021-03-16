/**
     (c) Copyright 2005
     Author: Ola Skavhaug
     Simula Research Laboratory AS
     
     This file is part of DegMatSparse.

     DegMatSparse is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.

     DegMatSparse is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with DegMatSparse; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

**/


#include "MapMatSparse.h"
#include <iostream>
#include <fstream>

MapMatSparse:: MapMatSparse(){}

double& MapMatSparse:: operator()(int i, int j) {
    pair<int,int> p(i,j) ;
    return values[p];
}

void MapMatSparse:: save(string& filename, string& name){
    ofstream oput (filename.c_str());
    map<pair<int, int>, double>::iterator iter;
    map<pair<int, int>, double>::iterator end = values.end();
    int i, j;
    int i_prev = -1;
    pair<int,int> p;
    end--;
    oput << "datatype = 'real';" << endl << endl;
    oput << "nrows = " << (*end).first.first+1 << "; ncolumns = " << (*end).first.second +1<< ";"<<endl;
    oput << "nentries= " << values.size() << ";" << endl << endl;
    oput <<  name << " = spalloc(nrows,ncolumns,nentries);"<<endl<<endl;
    oput << "% Data:" << endl << endl;

    end++;
    for (iter = values.begin(); iter !=end; iter++) {
        i = (*iter).first.first;
        j = (*iter).first.second;
        oput << name << "("<<i+1<<","<<j+1<<") = "<< (*iter).second<<";"<<endl; 
    }

}

void MapMatSparse:: load(char* filename) {
    int* i = new int(); 
    int* j = new int(); 
    float f[1];
    ifstream input (filename);
    char buffer[256];
    string content;
    input.getline(buffer, 255);
    while (!input.eof()){
        content = buffer;
        if (sscanf(content.c_str(), "A(%d,%d) = %f", i, j, f) >= 3) {
                (*this)((*i)-1,(*j)-1) = f[0];
        } else {
        }
        input.getline(buffer, 255);
    }


}

void MapMatSparse:: prod(valarray<double>* x, valarray<double>* b) {
    map<pair<int, int>, double>::iterator iter;
    map<pair<int, int>, double>::iterator end = values.end();
    int i, j;
    int i_prev = -1;
    pair<int,int> p;
    for (iter = values.begin(); iter !=end; iter++) {
        i = (*iter).first.first;
        j = (*iter).first.second;
        if (i != i_prev){ 
            i_prev = i;
            (*b)[i] = 0.0;
        }
        (*b)[i] += (*iter).second*(*x)[j];
    }
}
