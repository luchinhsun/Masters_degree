function x = tam(A,d)

    a(1) = 0.0;
    c(length(d)) = 0.0;
    for i = 1:length(d)-1
        c(i) = A(i,1+i);
        a(i+1) = A(i+1,i);
    end

    for i = 1:length(d)
        b(i) = A(i,i);
    end

    c_new(1) = c(1)/b(1);
    for i = 2:length(d)-1
        c_new(i) = c(i)/(b(i)-a(i)*c_new(i-1));
    end
    
    d_new(1) = d(1)/b(1);
    for i = 2:length(d)
        d_new(i) = (d(i)-a(i)*d_new(i-1))/(b(i)-a(i)*c_new(i-1));
    end

    x(length(d)) = d_new(length(d));
    for i = length(d)-1:-1:1
        x(i) = d_new(i)- c_new(i)*x(i+1);
    end
    
end