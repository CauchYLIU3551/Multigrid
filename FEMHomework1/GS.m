## Copyright (C) 2020 ChengyuLIU
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
## @deftypefn {} {@var{retval} =} GS (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: ChengyuLIU <ChengyuLIU@LAPTOP-732BUELF>
## Created: 2020-11-19

function x = GS (A, x, b, num)

n=size(A);

x0=x;
x0(1)=2*x(1);

while num>0||max(max(x-x0))>0.001
  x0=x;
  for i=1:n
    x(i)=(b(i)-A(i,:)*x+A(i,i)*x(i))/A(i,i); 
  end
  num=num-1;
end


end
