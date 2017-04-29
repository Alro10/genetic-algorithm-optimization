function elev=testfunction(loc)
% Cost function for the longitude-latitude example

[m,n]=size(loc);
for i=1:m
x=loc(i,1);
if x>2
x=2;
elseif x<-1
x=-1;
end
y=loc(i,2);
if y>2
y=2;
elseif y<-1
y=-1;
end
elev(i)=x*sin(4*pi*x)-y*sin(4*pi*y+pi)+1;
end