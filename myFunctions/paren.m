function y = paren(x, varargin)
% funtion handle to be able to index without assigning to variable
y = x(varargin{:});
end

