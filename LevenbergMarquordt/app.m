function app()
    syms x y;
    f1 = (x-9)^2+(y-6)^2+7;
    x0 = [9; 6];
    I =  eye(2);
    b0 = I*x0;
    deltaF = gradient(f1, [x y]);    
    
end

function f = problem(x)
    f = (x(1)-9)^2+(x(2)-6)^2+7;
end

function f1 = finit()
    syms x y;
    f1 = (x-9)^2+(y-6)^2+7;
end


