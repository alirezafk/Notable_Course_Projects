function Code=grayCodes(n)      %This function recieves input n and its output is a 2^n*n matrix which its ith row corresponds to the ith gray code with length n
    if(n==1)
        Code=[0;1];
    else
        Code=zeros(2^n,n);
        A=grayCodes(n-1);
        for i=1:2^(n-1)
            Code(i,:)=[0 A(i,:)];
        end
        for i=1:2^(n-1)
            Code(i+2^(n-1),:)=[1 A(2^(n-1)-i+1,:)];
        end
    end
end