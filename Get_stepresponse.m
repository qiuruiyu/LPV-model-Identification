function [StepR] = Get_stepresponse(w,span)
    K=0.6+w^3;
    T=3+0.5*w^4;
    num=K;
    den=[T,1];
    Ts=1;
    sysc=tf(num,den);
    sysd=c2d(sysc,Ts);
    StepR=step(sysc,1:span);
%     u = 0.5*ones(span, 1);
%     StepR = lsim(sysc, u, 1:span);
end
