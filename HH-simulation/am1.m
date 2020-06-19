function a=am1(v1) %Alpha for Variable m
    a=0.1*(v1+35)/(1-exp(-(v1+35)/10));
end