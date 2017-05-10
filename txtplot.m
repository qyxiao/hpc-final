
Npoints = 4;
xLength = Npoints; 
yLength = Npoints;

fin = fopen('outputComplex.txt', 'r');

A = fscanf(fin,'%f %f',[2*xLength*yLength 1]);

t = reshape(A,2,[]);

TheMatrix = reshape(complex(t(1,:),t(2,:)), xLength, yLength);

fclose(fin);