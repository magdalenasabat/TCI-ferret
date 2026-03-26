function [POSX,POSY] = myMakeCoords(X_vec, Y_vec, sz)



POSX = zeros(sz);
POSY = zeros(sz);

POSX(:,1) = linspace(X_vec(1),X_vec(2),size(POSX,1));
POSX(size(POSX,1),:) = linspace(X_vec(2),X_vec(3),size(POSX,2));
POSX(:,size(POSX,2)) = linspace(X_vec(4),X_vec(3),size(POSX,1));
POSX(1,:) = linspace(X_vec(1),X_vec(4),size(POSX,2));

for i=1:size(POSX,2)-2
    POSX(:,i+1) = linspace(POSX(1,i+1),POSX(end,i+1),size(POSX,1));
end

POSY(:,1) = linspace(Y_vec(1),Y_vec(2),size(POSY,1));
POSY(size(POSY,1),:) = linspace(Y_vec(2),Y_vec(3),size(POSY,2));
POSY(:,size(POSY,2)) = linspace(Y_vec(4),Y_vec(3),size(POSY,1));
POSY(1,:) = linspace(Y_vec(1),Y_vec(4),size(POSY,2));

for i=1:size(POSY,2)-2
    POSY(:,i+1) = linspace(POSY(1,i+1),POSY(end,i+1),size(POSY,1));
end

POSY = floor(flipud(POSY));
POSX = floor(flipud(POSX));
end