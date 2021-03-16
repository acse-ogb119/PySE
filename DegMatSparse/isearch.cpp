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

#include "isearch.h"
#include <valarray>

int interpolationSearch(int* vec, int key, int low, int high) {
  // Interpolation search for the index of "key"
    double low_diff, range_diff, count_diff;
    while ( vec[high] >= key && key > vec[low] ) {
        low_diff = (double)key - vec[low];
        range_diff = (double)vec[high] - vec[low];
        count_diff = (double)high - low;
        int range = (int)( low_diff / range_diff * count_diff + low );
        if ( key > vec[range] )
            low = range + 1;
        else if ( key < vec[range] )
            high = range - 1;
        else
            low = range;
    } 
    if ( key == vec[low] ) 
        return low;
    else 
        return -1;
};


int interpolationSearch(valarray<int>& vec, int key, int low, int high) {
  // Interpolation search for the index of "key"
    double low_diff, range_diff, count_diff;
    while ( vec[high] >= key && key > vec[low] ) {
        low_diff = (double)key - vec[low];
        range_diff = (double)vec[high] - vec[low];
        count_diff = (double)high - low;
        int range = (int)( low_diff / range_diff * count_diff + low );
        if ( key > vec[range] )
            low = range + 1;
        else if ( key < vec[range] )
            high = range - 1;
        else
            low = range;
    } 
    if ( key == vec[low] ) 
        return low;
    else 
        return -1;
};


