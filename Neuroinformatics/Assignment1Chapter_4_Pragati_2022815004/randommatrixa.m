function a = randommatrixa(nrows, ncols)
    a = rand(nrows, ncols);
    disp(a);
    
    for i = 1:nrows
        for j = 1:ncols
            k = a(i, j);
            
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
                fprintf('The %d %s row and %d %s column has a value of %.3f and is bigger than 0.5.\n', i, nrow_suffix, j, ncol_suffix, k);
            else
                fprintf('The %d %s row and %d %s column has a value of %.3f and is not bigger than 0.5.\n', i, nrow_suffix, j, ncol_suffix, k);
            end
        end
    end
end
