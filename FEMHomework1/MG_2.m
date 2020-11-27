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
## @deftypefn {} {@var{retval} =} MG_2 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: ChengyuLIU <ChengyuLIU@LAPTOP-732BUELF>
## Created: 2020-11-19

function x = MG_2 (A, b)
   n=size(A);
   I12=zeros((n+1)/2,n);
   for i=1:(n+1)/2-1
     I12(i,2*i)=1/2;
     I12(i,2*i-1)=1/4;
     I12(i,2*i+1)=1/4;
   end
   I12(end,end)=1;
   %I21=zeros(n,n/2);
   %for i=1:n/2-1
   %  I21(2*i,i)=1;
   %  I21(2*i-1,i)=1/2;
   %  I21(2*i+1,i)=1/2;
   %end  
   %I21(end,end)=1;
   I21=I12'*2;
   
   x=zeros(n,1);   
   x0=ones(n,1);

   while max(max(abs(x0-x)))>0.0000001
     x0=x;
     x=GS(A,x,b,10);
     r=A*x-b;
     tmp=(n+1)/2;
     e=zeros(tmp,1);
     e=GS(I12*A*I21,e,I12*r,20);
     x=x-I21*e;
   end
   
endfunction
