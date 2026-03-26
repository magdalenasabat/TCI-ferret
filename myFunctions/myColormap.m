function cmap = myColormap()
%%
load('my_cmap.mat')

beg = original(125:end,:);
en = original(1:50,:);

necmap = [beg; en];

test = 1;

if test

   M = repmat([1:1:256],256,1);
   figure();
   imagesc(M);
   colormap(necmap);
   colorbar
    
    
end
%%
my_cmap = necmap;
end