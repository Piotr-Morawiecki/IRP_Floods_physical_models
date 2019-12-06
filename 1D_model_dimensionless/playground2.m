n = 1.282;
m = 1 - 1 / n;
s = (0:100)/10000;
func = (1-s.^(n-1).*(1+s.^n).^(-m)).^2./(1+s.^n).^(m./2);
approx = (1-s.^(n-1)+2.*m.*s.^(2.*n-1)+s.^(2.*n-2)-n.*s.^(3.*n-2)+m.*s.^(2.*n-1)-m./2.*s.^(3.*n-2)); 

plot(s, func);
hold on
plot(s, approx);
hold off