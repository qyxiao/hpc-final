
Npoints = 1024;
filesNum = 16;
xLength = Npoints; 
yLength = Npoints;


for i = 0:(filesNum-1)
   
    fin = fopen(['outputMPI',num2str(i),'.txt'], 'r');
    A = fscanf(fin,'%f %f',[2*(xLength/filesNum)*yLength 1]);
    t = reshape(A,2,[]);
    TheMatrix = reshape(complex(t(1,:),t(2,:)), 1, (xLength/filesNum)*yLength);
    if(i==0)
        output = TheMatrix;
    else
        output = [output TheMatrix];
    end
    
    fclose(fin);
end

output = reshape(output,xLength,yLength);
figure(1);clf
pcolor(real(output));shading flat;colorbar 


