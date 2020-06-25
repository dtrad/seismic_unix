t=-128:127;
figure(1)
plot(t,t.^2,'o',t,abs(t),t,log(1000+t.^2))
