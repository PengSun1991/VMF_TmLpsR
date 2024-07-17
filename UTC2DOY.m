function DOY=UTC2DOY(y,m,d)

JD0=fix(365.25*(y-1))+fix(30.6001*(13+1))+1+1720981.5;

if m<=2
    y=y-1;m=m+12;
end

JD1=fix(365.25*y)+fix(30.6001*(m+1))+d+1720981.5;
DOY=JD1-JD0+1;
end

