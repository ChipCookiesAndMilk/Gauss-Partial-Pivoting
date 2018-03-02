% Matrix definition
% a = [a11,a12,a13,a14; a21,a22,a23,a24; a31,a32,a33,a34; a41,a42,a43,a44];
function GaussPartialPivoting(a) 
    %a(r,c) = original(r,c);
    % r - rows, 	y
    % c - colums, 	x
    [r, c] = size(a);
    
    if c-1 == r
        % VARIABLES
            % fill up the array x with zeros
            x       = zeros(r,1); 
            % first pivot
            value 	= a(1,1);
            % if a better pivot is found
            flag 	= -1;
            % position in x array
            n       = r;
            % auxiliar variable
            sub 	= 0;
        % GAUSS ALGORITHM
        fprintf("\nInserted matrix\n");
        display(a);
        for j=1:r-1
            %% PART A, pivoting if necessary
            for z=j:r
                if value<a(z,j) 
                    flag = z;
                    value = a(z,j);
                end
            end
            if flag>-1
                t           = a(j,:);
                a(j,:)      = a(flag,:);
                a(flag,:)   = t;
            end
            %% Part B, Elimination pass on all rows i>k
            for i = j+1:r
                a(i,:) = a(i,:)-a(j,:)*(a(i,j)/a(j,j));
            end
            flag = -1;
        end
        fprintf("\nNew matrix\n");
        disp(a);
        disp(" ");
        %% Part C, Blackwards Substitution  
        % get first x
        x(n) = a(r,c)/a(r,c-1);
        fprintf("Resulting x\n\tx(%d): %d\n",n,x(n));

        % for the rest of the x's positions
        for i=n-1:-1:1
            % 1.1 Substraction
            for j=c-1:-1:i+1
                sub = sub - (a(i,j)*(x(j)));
            end
            % 1.2 Obtain x
            x(i) = (sub + a(i,c))/a(i,i);
            sub = 0;
            fprintf("\tx(%d): %d\n",i,x(i));
        end
    else
        fprintf("\n\tThe matriz must be extended, the last column should contain\nthe [B] values\n");
    end
end