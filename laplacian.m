function arrOut = laplacian(arrIn)
% Calculates laplacian for a matrix using a 3x3 convolution with edge wrapping
arrOut = -1*arrIn + ...
    0.2*(circshift(arrIn,1)+circshift(arrIn,-1)+circshift(arrIn,-1,2)+circshift(arrIn,1,2)) + ...
    0.05*(circshift(arrIn,[1 1])+circshift(arrIn,[1 -1])+circshift(arrIn,[-1 1])+circshift(arrIn,[-1 -1]));
end
