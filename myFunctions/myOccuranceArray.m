function output = myOccuranceArray(A)
    output = zeros(size(A));

    unique_elements = unique(A);
    for i = 1:length(unique_elements)
        element = unique_elements(i);
        indices = (A == element);
        counts = sum(indices);
        output(indices) = 1:counts;
    end
end
