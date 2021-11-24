function [B, e]=sigmoidFit(ii,x,y)
x = x(ii); y = y(ii); 
options = optimset('MaxFunEvals',10000);
[B, e] = fminsearch(@(z) mysigmoidfit(z,x,y),[1.0;0.2; min(y); max(y)],options);
% if abs(B(3))>abs(B(4))
%     disp("sigmoid fits other way around")
%     B(3) = (B(3) + B(4))/2;
%     B(4) = (B(3) + B(4))/2;
%     e = mysigmoidfit(B,x,y,ve)
% end
return
end