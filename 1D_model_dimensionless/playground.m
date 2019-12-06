h = ((1:10000)-5000)/5000;
hold on 
plot(h, computeKrDerivative (h, mod));
hold off
plot(h, (computeKr (h+0.0001, mod)-computeKr (h, mod))/0.0001);



[z,y] = ode45(@fun,[0 20],[2; 0]);

plot(t,y(:,1),'-o',t,y(:,2),'-o')
title('Solution of van der Pol Equation (\mu = 1) with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')

function dy = fun(z,y)
    % y1 = h
    % y2 = q
    % dh/dz=
    dy = [y(2); (1-y(1)^2)*y(2)-y(1)];
end