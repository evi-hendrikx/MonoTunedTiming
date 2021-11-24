function [b_lower, b_upper, xfit, yfit, xfitR2, yfitR2] = monoTuned_make_sigmoid_CIs_anonymous(x_steps,y_values)

[B] = bootstrp(10000,@(x) sigmoidFit(x, x_steps, ...
    y_values), 1:length(x_steps));

% 1000 steps for plots
xfit = linspace(min(x_steps),max(x_steps),1000)';
for n=1:size(B, 1)
    yfitAll(n,:) = B(n,4)-normcdf2(xfit,B(n,1),B(n,2)).*(B(n,4)-B(n,3));
end

% calculate 95% CI of evaluated parameter fits
pct1 = 100*0.05/2;
pct2 = 100-pct1;
b_lower=prctile(yfitAll,pct1);
b_upper=prctile(yfitAll,pct2);

%figure; plot(xfit, b_lower, 'r'); hold on; plot(xfit, b_upper, 'g'); plot(xfit, yfit, 'k'); %Use these values, but no need to plot them yet

% make sigmoid
main_sigmoid_fit = sigmoidFit(1:length(x_steps), x_steps, y_values);

% for plot
yfit = main_sigmoid_fit(4)-normcdf2(xfit,main_sigmoid_fit(1),main_sigmoid_fit(2)).*(main_sigmoid_fit(4)-main_sigmoid_fit(3));

% for calculating the residuals
xfitR2 = x_steps;
yfitR2 = main_sigmoid_fit(4)-normcdf2(xfitR2,main_sigmoid_fit(1),main_sigmoid_fit(2)).*(main_sigmoid_fit(4)-main_sigmoid_fit(3));

end