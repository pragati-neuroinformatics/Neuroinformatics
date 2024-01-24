%% Question 10 (generating the random matrix 32*3 when (1 for bigger than 0.5; 0 for smaller than 0.5)
rows = 32;
columns = 3;

resultMatrix = randommatrixC(rows,columns); % performed using randommatrixC function 

for i = 1:rows
    for j = 1:columns
        if resultMatrix(i, j) > 0.5
            resultMatrix(i, j) = 1;
        else
            resultMatrix(i, j) = 0;
        end
    end
end
fprintf('Result Matrix (Row, Column, Result):\n');
disp(resultMatrix);

%% Question 11 (Assigned the appropriate variable labels and made it readble in spreadsheet version)
disp (resultMatrix)

variableLabels = {'Row_Index', 'Coloumn_Index', 'Test_Results'};

% Combining variable labels and data into a table
dataTable = array2table(resultMatrix, 'VariableNames', variableLabels);

% Saving the table to a tab-delimited file(TAB)
writetable(dataTable, 'RandomMatrixScriptC.txt', 'Delimiter', '\t');

