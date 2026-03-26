function [POSX, POSY, imageData,distances ] = myCalculatePosition(indiv, sessN)

sz = [4,9];

if strcmpi(indiv, 'murols') ;
    X = [  246 191 303 348];
    Y = [  278 199 120 175];
%     X = [  196 260 256 199];
%     Y = [  117 117 290 290];
    t = Tiff('brain_L_new.tiff','r');
    imageData = read(t);
    [POSX,POSY] = myMakeCoords_v2(X, Y, sz); % calculate for each electrode 
    dist = load('ferret_brain_anatomy_distancefromcenter_L.mat');
    dist = dist.d_L;
elseif strcmpi(indiv, 'frinault');
    if all( sessN <= 17 );
        X = [92 117 252  240 ];
        Y = [ 317 251 279  356] ;
        t = Tiff('brain_R_new.tiff','r');
        imageData = read(t);
        [POSX,POSY] = myMakeCoords_v2(X, Y, sz); % calculate for each electrode 
        dist = load('ferret_brain_anatomy_distancefromcenter_R.mat');
        dist = dist.d_R;
    elseif  all(sessN > 17 &  sessN < 23);
        X = [213 233 80 75];
        Y = [222 285 325 250 ];
        t = Tiff('brain_L_new.tiff','r');
        imageData = read(t);
        [POSX,POSY] = myMakeCoords_v2(X, Y, sz); % calculate for each electrode 
        dist = load('ferret_brain_anatomy_distancefromcenter_L.mat');
        dist = dist.d_L;
    else strcmpi(indiv, 'frinault'); %right plus left
        X1 = [ 92 117 252  240] + 383;
        Y1 = [ 317 251 279  356] ;
        X = [213 233 80 75];
        Y = [222 285 325 250 ] ;

        t1 = Tiff('brain_L_new.tiff','r');
        imageData1 = read(t1);
        t2 = Tiff('brain_R_new.tiff','r');
        imageData2 = read(t2);

        imageData = [imageData1 imageData2];
        [POSX1,POSY1] = myMakeCoords_v2(X, Y, sz); % calculate for each electrode 
        [POSX2,POSY2] = myMakeCoords_v2(X1, Y1,sz); % calculate for each electrode 

        POSX = [POSX1 POSX2];
        POSY = [POSY1 POSY2];
        d_L = load('ferret_brain_anatomy_distancefromcenter_L.mat');
        d_R = load('ferret_brain_anatomy_distancefromcenter_R.mat');
        dist = [d_L.d_L d_R.d_R]; 
    end
elseif strcmpi(indiv, 'manigodine')   
%     X = [86 100 112 99];
%     Y = [120 115 144 148 ];
    X = [7 21 165 146 ];
    Y = [194 146 202 247  ];
    t = Tiff('brain_R_new.tiff','r');
    imageData = read(t);
    [POSX,POSY] = myMakeCoords_v2(X, Y, sz); % calculate for each electrode 
    dist = load('ferret_brain_anatomy_distancefromcenter_R.mat');
    dist = dist.d_R;
end

distances_tmp = [];
for i =1: length(POSX(:)); distances_tmp(i) = dist(POSY(i),POSX(i)); end
distances = reshape(distances_tmp, size(POSX));
            
end
