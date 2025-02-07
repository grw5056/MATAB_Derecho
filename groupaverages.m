function [avg]=groupaverages(input_vector,depth)
    m=1;
    input_vector=sort(input_vector);
    while any(input_vector~=inf)
        if input_vector(1)==inf
            input_vector=cat(2,input_vector(2),input_vector(3:end),input_vector(1));
        else
            for n=2:length(input_vector)
                avg_num=1;
                avg(m)=input_vector(1);
                if abs(input_vector(1)-input_vector(n))<=depth/1000
                    avg(m)=avg(m)+input_vector(n);
                    avg_num=avg_num+1;
                    input_vector(n)=inf;
                end
            end
            input_vector(1)=inf;
            avg(m)=avg(m)/avg_num;
            m=m+1;
        end
    end
