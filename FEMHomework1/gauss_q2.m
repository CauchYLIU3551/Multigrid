## Copyright (C) 2020 63002
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} gauss_q2 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: 63002 <63002@DESKTOP-72PSU6T>
## Created: 2020-11-27

function F = gauss_q2 (a, b)

r1=0.7745966692;
r2=0;
r3=-r1;
c1=5/9;
c2=8/9;
c3=5/9;
x=[r1;r2;r3];
c=[c1,c2,c3];
%b-a/2 * f(((b-a)t+(b+a))/2);
x=((b-a).*x+(b+a))./2;
%F=(b-a)/2*(x.^6-x.^2.*sin(2*x));
F=(b-a)/2*(-exp(x).*(cos(x)-2*sin(x)-x.*cos(x)-x.*sin(x)).*(b-x)/(b-a));
F=c*F;
endfunction
