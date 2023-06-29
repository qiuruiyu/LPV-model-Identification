function [StepR] = Get_stepresponse(w,span)
    K = 0.3+w^2;
    T = 3+0.5*w^4;
    num=K;
%     den=[T,0.1];
    den = [1,2,1];
    Ts=1;
    sysc=tf(num,den);
    sysd=c2d(sysc,Ts);
    StepR=step(sysc,1:span);
%     u = 0.5*ones(span, 1);
%     StepR = lsim(sysc, u, 1:span);
end
