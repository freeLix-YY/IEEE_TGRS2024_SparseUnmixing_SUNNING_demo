function [A] = BuildLibraryCuprite()
    load USGS_1995_Library.mat;

    %  order bands by increasing wavelength
    [dummy index] = sort(datalib(:,1));
    A = datalib(index,4:end);
    names = names(4:end,:);

    % order the columns of A by increasing angles 
    [A, index, angles] = sort_library_by_angle(A);
    names = names(index',:);

end