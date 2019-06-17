clear all, clc

X = load('Binary.txt');
size(X)

T = round(X*X');
size(T)
FT = fopen('Tree.txt', 'w')
for i=1:length(T)
    for j=1:length(T)
        fprintf(FT, '%i ', T(i, j));
    end
    fprintf(FT, '\n');
end
%save('Tree.txt', 'T', '-ascii')

F = round(X'*X);
size(F)
FF = fopen('Fungi.txt', 'w')
for i=1:length(F)
    for j=1:length(F)
        fprintf(FF, '%i ', F(i, j));
    end
    fprintf(FF, '\n');
end
%save('Fungi.txt', 'F', '-ascii')
fclose all