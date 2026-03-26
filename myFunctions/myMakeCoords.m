function [POSX,POSY] = myMakeCoords(theta,mov_X,move_Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates coordinates of the electrodes given what we know about the
% array and given the user defined coordinates. It achors to X1 and Y1 and
% uses the rest to determine the directions
%
% Magdalena Sabat, 15/09/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% electrode pitch is 400 micrometers -> 20 pixels
POSX_orig = [30,50,70,90,110,130,150,170,190; ...
    20,40,60,80,100,120,140,160,180;...
    30,50,70,90,110,130,150,170,190;...
    20,40,60,80,100,120,140,160,180] ;

% therefore the row height is ~346 micrometers -> ~340 -> 17 pixels
POSY_orig  = [1,1,1,1,1,1,1,1,1;...
    18,18,18,18,18,18,18,18,18;...
    35,35,35,35,35,35,35,35,35;...
    52,52,52,52,52,52,52,52,52] ;

mask = logical( [1   0   0   0   0   0   0   0   0 ; ....
   0   0   0   0   0   0   0   0   1 ; ...
   1   0   0   0   0   0   0   0   0 ; ...
   0   0   0   0   0   0   0   0   1]);



% Create rotation matrix
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

% Rotate all points
tmp = [POSX_orig(:),POSY_orig(:)]';

POS = R*tmp;
POSX= reshape(floor(POS(1,:)+ mov_X),size(POSY_orig));
POSY = reshape(floor(POS(2,:)+ move_Y),size(POSX_orig));

% remove the masked values
% POSX(mask(:)) = 0;
% POSY(mask) = NaN;
end
