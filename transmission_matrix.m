function [A] = transmission_matrix(Phi_s, Alpha, Lx_p, Ly_p, Ly_w, Np, XC, YC, Area, Ventilation)

% PEN-TO-PEN TRANSMISSION MATRIX
% Authors: Maryam Safari, Christian Fleming, Jason A. Galvis,
% Aniruddha Deka, Felipe Sanchez, Gustavo Machado & Chi-An Yeh
% Email: msafari2@ncsu.edu
% North Carolina State University, Raleigh, North Carolina, USA.
% Department of Mechanical and Aerospace Engineering
% Department of Population Health and Pathobiology

% For details, please see:
%
% Ref:  M. Safari, C. Fleming, J. A. Galvis, A. Deka, F. Sanchez,
%       G. Machado, C.-A. Yeh,
%       "A CFD-informed barn-level livestock/swine disease dissemination
%       model and its use for ventilation optimization"
%       Epidemics Journal, 2025.

% Inputs:
% -- Phis_s       : Generation rate of the source
% -- Alpha        : Diffusion of the pathogen
% -- Lx_p         : Length of each pen
% -- Ly_p         : Width of each pen
% -- Ly_w         : Width of the walkway
% -- Np           : Number of pens per each row in the barn
% -- XC           : X coordinate of the cells center (Matrix)
% -- YC           : Y coordinate of the cells center (Matrix)
% -- Area         : Area of the cells (Matrix)
% -- Ventilation  : Ventilation setting as the BC in CFD solver

if any(size(XC)==1) || any(size(YC)==1)
    error('-- Cell coordinates must be a matrix')
end

if any(size(Area)==1)
    error('-- Area of the cells must be a matrix')
end
%% STEP 1 -- Generate indexing matrix for each cell

nx = size(XC, 1);
ny = size(YC, 2);

icv = nan(nx, ny);

count = 1 ;
for i = 1:nx
    for j = 1:ny
        icv(i, j) = count;
        count = count + 1;
    end
end

%% STEP 2 -- Define pen area according to generated index matrix

Pen_ID = cell(Np, 1) ;

for ipen = 1:2*Np

    if ipen <= Np
        Pen_X = [0, Lx_p, Lx_p, 0] + (ipen - 1)*Lx_p;
        Pen_Y = [0, 0, Ly_p, Ly_p];
    else
        Pen_X = [0, Lx_p, Lx_p, 0] + (ipen - Np - 1)*Lx_p;
        Pen_Y = [Ly_p + Ly_w, Ly_p + Ly_w, 2*Ly_p + Ly_w, 2*Ly_p + Ly_w];
    end

    ID = false(nx*ny, 1) ;
    for i = 1:nx
        for j = 1:ny
            if inpolygon(XC(i, j), YC(i, j), Pen_X, Pen_Y)
                ID(icv(i, j)) = true;
            end
        end
    end

    Pen_ID{ipen} = find(ID);
end

%% Step 3 -- Seed pathogen in each pen to construct the transmission matrix 
A = zeros(2*Np, 2*Np);

for isrc = 1:2*Np

    % Seed source in each pen, and form a source matrix
    Src = zeros(nx*ny, 1);

    area = reshape(Area', nx*ny, 1);
    Src(Pen_ID{isrc}) = Phi_s*area(Pen_ID{isrc});
    Source = Src(icv);

    % Run the CFD simulation to obtain the concentration field of airborne
    % pathogens with the prescribed source

    [Phi] = CFD_solver(Source, Alpha, Ventilation, XC, YC);

    % Render one column of the transmission matrix via each CFD simulation
    phi =  reshape(Phi',  nx*ny, 1) ;

    for ipen = 1:2*Np
        A(ipen, isrc) = sum(phi(Pen_ID{ipen}))/(Lx_p*Ly_p) ;
    end

end

end % function

