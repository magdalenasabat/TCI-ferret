function myTable2Latex(T,savedir)% Sample table
% chatGPT generated Table "T" to Latex code generation to a path in savedir
% magdalena.sabat@ens.psl.eu
arguments
    T table
    savedir char
end

    % Open a file to write LaTeX code
    fid = fopen(savedir, 'w');
    fprintf(fid, '\\begin{table}[h]\n\\centering\n\\begin{tabular}{|c|c|}\n');
    fprintf(fid, '\\hline\n');

    % Write the column headers
    fprintf(fid, '%s & %s \\\\ \\hline\n', T.Properties.VariableNames{1}, T.Properties.VariableNames{2});

    % Write the data rows
    for i = 1:height(T)
        fprintf(fid, '%d & %.2f \\\\ \\hline\n', T{i,1}, T{i,2});
    end

    % Close the LaTeX table
    fprintf(fid, '\\end{tabular}\n\\caption{Sample Table}\n\\end{table}\n');
    fclose(fid);
end