function e=mysigmoidfit(z,x,y)
e=sum((y-(z(4)-(normcdf2(x,z(1),z(2)).*(z(4)-z(3))))).^2);
return
end