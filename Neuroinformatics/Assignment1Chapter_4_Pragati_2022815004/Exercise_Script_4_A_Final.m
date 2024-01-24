% Exercise for Script A 
% Question_1
nrows=4;
ncols=8;
Rmatrix = rand(4,8); 
disp (Rmatrix)
%%
% Question_2 and Question_3 
    a = rand(nrows,ncols);
    disp (a);
    for i = 1:nrows
        for j=1:ncols
            if a(i,j)>0.5
            fprintf ('The %dth row and %dth column has a value of %.3f and is bigger than 0.5.\n', i, j, a(i, j));
            else
            fprintf ('The %dth row and %dth column has a value of %.3f and is not bigger than 0.5.\n', i, j, a(i, j));
            end
        end
    end
     %% Question 4 (to add exceptions to print out 1st, 2nd and 3rd). 
       nrows=4;
       ncols=8;
        a = rand(nrows,ncols);
        disp (a);
        for i = 1:nrows
            for j=1:ncols
                k= a(i,j);
               
                % changing_row_suffix
                nrow_suffix = 'th';
                if i == 1
                   nrow_suffix = 'st';
                elseif i == 2
                    nrow_suffix = 'nd';
                elseif i == 3
                    nrow_suffix = 'rd';
                end
                % changing_column_suffix
                ncol_suffix = 'th';
                if j == 1
                   ncol_suffix = 'st';
                elseif j == 2
                    ncol_suffix = 'nd';
                elseif j == 3
                    ncol_suffix = 'rd';
                end
                   
                if k > 0.5
                fprintf ('The %d %s row and %d %s column has a value of %.3f and is bigger than 0.5.\n', i,nrow_suffix,j, ncol_suffix,k);
                else
                fprintf ('The %d %s row and %d %s column has a value of %.3fand is not bigger than 0.5.\n', i,nrow_suffix, j, ncol_suffix,k);
                end
            end
        end
             

     %% Question_5 (created a function randommatrixa to call from the command line when 
     % we have 2 inputs: any number of rows and number of coloumns
a = randommatrixa(4,8);
