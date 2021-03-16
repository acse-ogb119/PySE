%
%     (c) Copyright 2005
%     Author: Ola Skavhaug
%     Simula Research Laboratory AS
%     
%     This file is part of DegMatSparse.
%
%     DegMatSparse is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
%
%     DegMatSparse is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with DegMatSparse; if not, write to the Free Software
%     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%

datatype = 'real';

nrows = 6; ncolumns = 6;
nentries= 7;

X = spalloc(nrows,ncolumns,nentries);

% Data:

A(1,1) = 0
A(1,4) = 19
A(2,2) = 1
A(3,3) = 2
A(4,4) = 3
A(5,5) = 4
A(6,6) = 5
