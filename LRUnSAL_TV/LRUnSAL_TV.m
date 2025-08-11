% This is the code re-implementation of LRUnSAL-TV
% For more information, see the following:
% Website: https://ieeexplore.ieee.org/document/8767806

function [X,res,rmse] = LRUnSAL_TV(M,Y,varargin)

%% --------------- Description ---------------------------------------------
%
%  LRUnSAL-TV solves the following nuclear norm + TV + l_1-norm optimization problem:
%
%     Definitions:
%
%      A  -> L * n; Mixing matrix (Library)
%      X  -> n * N; collection of N fractional vectors; each column  of X
%                   contains the fraction of a correspondent  pixel
%
%      Optimization problem:
%
%    min  (1/2) ||Y - AX - S||^2_F  + lambda_1 ||X||_{*}
%    X,S                         + lambda_tv ||LX||_{1,1} + lambda_2 ||S||_1;
%
%
%    where
%
%        (1/2) ||Y - AX - S||^2_F is a quadratic data misfit term
%
%        ||X||_{*} = nuclear norm
%
%        ||LX||_{1,1} is the TV (non-isotropic or isotropic regularizer)
%
%        ||S||_1 is the sparse noise term
%
%
%         L is a linear operator that computes the horizontal and the
%         vertical differences on each band of X.  Let Lh: R^{n*N}-> R^{n*N}
%         be a linear operator that computes the horizontal first order
%         differences per band. LhX computes a matrix of the same size of X
%         (we are assuming cyclic boundary), where [LhX](i,j) = X(i,h(j))-X(i,j),
%         where h(j) is the index of pixel on the right hand side of j.
%
%         For the vertical differnces, we have a similar action of Lv:
%         [LvX](i,j) = X(v(i),j)-X(i,j), where  v(i) is the index of pixel
%         on the top hand side of j.
%
%         We consider tow types of Total variation:
%
%         a) Non-isotropic:  ||LX||_{1,1} := ||[Lh; Lv]X||_{1,1}
%
%         b) Isotropic:  ||LX||_{1,1}  := ||(LhX, LvX)||_11,
%             where   |||(A,B)||_{1,1} := |||sqrt(A.^2 + B.^2)||_{1,1}
%
%
% -------------------------------------------------------------------------
%
%
%
%    CONSTRAINTS ACCEPTED:
%
%    1) Positivity X(:,i) >= 0, for i=1,...,N
%    2) Sum-To-One sum( X(:,i)) = 1, for for i=1,...,N
%
%
%    NOTES:
%
%       1) If X is a matrix and lambda_TV = 0, SUNSAL_TV solves
%           columnwise independent optimizations.
%
%       2) If both the Positivity and Sum-To-One constraints are active,
%          then we have ||X||_{1,1} = n and therefore this regularizer
%          is useless.
%
%
%% -------------------- Line of Attack  -----------------------------------
%
%  LRUnSAL_TV solves the above optimization problem by introducing a variable
%  splitting and then solving the resulting constrained optimization with
%  the augmented Lagrangian method.
%
%
%   The initial problem is converted into
%
%    min  (1/2) ||Y - AX - S||^2_F  + i_R_+(X)
%    X,S                            + i_S(X)
%                                   + lambda_1  ||X||_{*}
%                                   + lambda_tv ||LX||_{1,1} 
%                                   + lambda_2 ||S||_{1};
%
%
%   where i_R_+ and i_S are the indicator functions of the set R_+ and
%   the probability simplex, respecively, applied to the columns ox X.
%
%
%  Then, we apply the following variable splitting
%
%
%    min  (1/2) ||Y - V1 - S||^2     + i_R_+(V2)
%  X,V1, .... V7                     + i_S(V3)
%                                    + lambda_1  ||V4||_{*}
%                                    + lambda_tv ||V6||_{1,1} 
%                                    + lambda_2 ||S||_{1};
%
%     subject to:  AX   = V1
%                  X    = V2
%                  X    = V3
%                  X    = V4
%                  X    = V5
%                  HV5  = V6
%
%
%
%
% ------------------------------------------------------------------------
%%  ===== Required inputs =============
%
%  M - [L(observations) * n (variables)] system matrix (usually a library)
%
%  Y - matrix with  L(observation) x N(pixels).
%
%
%%  ====================== Optional inputs =============================
%
%
%  'LAMBDA_1' - regularization parameter for nuclear norm.
%               Default: 0;
%
%  'LAMBDA_TV' - regularization parameter for TV norm.
%                Default: 0;
%
%  'LAMBDA_2' - regularization parameter for l1 norm.
%                Default: 0;
%
%  'TV_TYPE'   - {'iso','niso'} type of total variation:  'iso' ==
%                isotropic; 'n-iso' == non-isotropic; Default: 'niso'
%
%  'IM_SIZE'   - [nlins, ncols]   number of lines and rows of the
%                spectral cube. These parameters are mandatory when
%                'LAMBDA_TV' is  passed.
%                Note:  n_lin*n_col = N
%
%
%  'AL_ITERS' - (double):   Minimum number of augmented Lagrangian iterations
%                           Default 100;
%
%
%  'MU' - (double):   augmented Lagrangian weight
%                           Default 0.001;
%
%
%
%  'POSITIVITY'  = {'yes', 'no'}; Default 'no'
%                  Enforces the positivity constraint: x >= 0
%
%  'ADDONE'  = {'yes', 'no'}; Default 'no'
%               Enforces the positivity constraint: x >= 0
%
%  'TRUE_X'  - [n (variables), N (pixels)] original data in matrix format.
%              If  the XT (the TRUE X) is inputted, then the RMSE is
%              ||X-XT||computed along the iterations
%
%
%  'VERBOSE'   = {'yes', 'no'}; Default 'no'
%
%                 'no' - work silently
%                 'yes' - display warnings
%
%%  =========================== Outputs ==================================
%
% X  =  [nxN] estimated  X matrix
%
%
% ----------------------------------------------------------------------

%%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 2
    error('Wrong number of required parameters');
end
% mixing matrix size
[LM,n] = size(M);
% data set size
[L,N] = size(Y);
if (LM ~= L)
    error('mixing matrix M and data set y are inconsistent');
end




%%
%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
%


% 'LAMBDA_1'
%  nuclear regularization
reg_nuclear = 0; % absent

% 'LAMBDA_TV'
%  TV regularization
reg_TV = 0; % absent
im_size = []; % image size
tv_type = 'niso'; % non-isotropic TV

% 'LAMBDA_2'
%  l1 regularization
reg_l1 = 0; % absent

% 'AL:ITERS'
% maximum number of AL iteration
AL_iters = 1000;

% 'MU'
% AL weight
mu = 0.001;

% 'VERBOSE'
% display only warnings
verbose = 'off';

% 'POSITIVITY'
% Positivity constraint
positivity = 'no';
reg_pos = 0; % absent

% 'ADDONE'
%  Sum-to-one constraint
addone = 'no';
reg_add = 0; % absent

%

% initialization
X0 = 0;

% true X
true_x = 0;
rmse = 0;

% tolerance for the primal and dual residues
tol = 1e-6;

%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------


%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'LAMBDA_1'
                lambda_1 = varargin{i+1};
                if lambda_1 < 0
                    error('lambda must be positive');
                elseif lambda_1 > 0
                    reg_nuclear = 1;
                end
            case 'LAMBDA_TV'
                lambda_TV = varargin{i+1};
                if lambda_TV < 0
                    error('lambda must be non-negative');
                elseif lambda_TV > 0
                    reg_TV = 1;
                end
            case 'LAMBDA_2'
                lambda_2 = varargin{i+1};
                if lambda_2 < 0
                    error('lambda must be positive');
                elseif lambda_2 > 0
                    reg_l1 = 1;
                end    
            case 'TV_TYPE'
                tv_type = varargin{i+1};
                if ~(strcmp(tv_type,'iso') | strcmp(tv_type,'niso'))
                    error('wrong TV_TYPE');
                end
            case 'IM_SIZE'
                im_size = varargin{i+1};
            case 'AL_ITERS'
                AL_iters = round(varargin{i+1});
                if (AL_iters <= 0 )
                    error('AL_iters must a positive integer');
                end
            case 'POSITIVITY'
                positivity = varargin{i+1};
                if strcmp(positivity,'yes')
                    reg_pos = 1;
                end
            case 'ADDONE'
                addone = varargin{i+1};
                if strcmp(addone,'yes')
                    reg_add = 1;
                end
            case 'MU'
                mu = varargin{i+1};
                if mu <= 0
                    error('mu must be positive');
                end
            case 'TOL'
                tol = varargin{i+1};    
            case 'VERBOSE'
                verbose = varargin{i+1};
            case 'X0'
                X0 = varargin{i+1};
            case 'TRUE_X'
                XT = varargin{i+1};
                true_x = 1;
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

% test for true data size correctness
if true_x
    [nr nc] = size(XT);
    if (nr ~= n) | (nc ~= N)
        error('wrong image size')
    end
end


% test for image size correctness
if reg_TV > 0
    if N ~= prod(im_size)
        error('wrong image size')
    end
    n_lin = im_size(1);
    n_col = im_size(2);
    
    % build handlers and necessary stuff
    % horizontal difference operators
    FDh = zeros(im_size);
    FDh(1,1) = -1;
    FDh(1,end) = 1;
    FDh = fft2(FDh);
    FDhH = conj(FDh);
    
    % vertical difference operator
    FDv = zeros(im_size);
    FDv(1,1) = -1;
    FDv(end,1) = 1;
    FDv = fft2(FDv);
    FDvH = conj(FDv);
    
    IL = 1./( FDhH.* FDh + FDvH.* FDv + 1);
    
    Dh = @(x) real(ifft2(fft2(x).*FDh));
    DhH = @(x) real(ifft2(fft2(x).*FDhH));
    
    Dv = @(x) real(ifft2(fft2(x).*FDv));
    DvH = @(x) real(ifft2(fft2(x).*FDvH));
    
end




%%
%---------------------------------------------
%  Constants and initializations
%---------------------------------------------

% number of regularizers
n_reg =  reg_nuclear + reg_pos + reg_add + reg_TV;

IF = inv(M'*M + n_reg*eye(n));

%%
%---------------------------------------------
%  Initializations
%---------------------------------------------

% no intial solution supplied
if X0 == 0
    X = IF*M'*Y;
end

% what regularizers ?
%  1 - data term
%  2 - positivity
%  3 - addone
%  4 - nuclear
%  5 - TV
%  6 - l1 (noise)


index = 1;

% initialize V variables
V = cell(1 + n_reg,1);

% initialize D variables (scaled Lagrange Multipliers)
D = cell(1 + n_reg,1);


%  data term (always present)
reg(1) = 1;             % regularizers
V{index} = M*X;         % V1
D{1} = zeros(size(Y));  % Lagrange multipliers

% next V
index = index + 1;
% POSITIVITY
if reg_pos == 1
    reg(index) = 2;
    V{index} = X;
    D{index} = zeros(size(X));
    index = index +1;
end
% ADDONE
if reg_add == 1
    reg(index) = 3;
    V{index} = X;
    D{index} = zeros(size(X));
    index = index +1;
end
% nuclear norm (SVT)
if reg_nuclear == 1
    reg(index) = 4;
    V{index} = X;
    D{index} = zeros(size(X));
    index = index +1;
end
% TV
% NOTE: V5, V6, D5, and D6 are represented as image planes
if reg_TV == 1
    reg(index) = 5;
    % V5
    V{index} = X;
    D{index} = zeros(size(X));
    
    % convert X into a cube
    U_im = reshape(X',im_size(1), im_size(2),n);
    
    % V6 create two images per band (horizontal and vertical differences)
    V{index+1} = cell(n,2);
    D{index+1} = cell(n,2);
    for i=1:n
        % build V6 image planes
        V{index+1}{i}{1} = Dh(U_im(:,:,i));   % horizontal differences
        V{index+1}{i}{2} = Dv(U_im(:,:,i));   % horizontal differences
        % build d7 image planes
        D{index+1}{i}{1} = zeros(im_size);   % horizontal differences
        D{index+1}{i}{2} = zeros(im_size);   % horizontal differences
    end
    clear U_im;
end


%%
%---------------------------------------------
%  AL iterations - main body
%---------------------------------------------
tol1 = sqrt((3*n+L)*N)*tol;
tol2 = sqrt((3*n+L)*N)*tol;

i=1;
res_p = inf;
res_d = inf;
while (i <= AL_iters) && ((abs(res_p) > tol1) || (abs(res_d) > tol2))
    
    if mod(i,10) == 1
        for j = 1:(n_reg+1)
            V10{j} = V{j};
        end
    end
    
    % solve the quadratic step (all terms depending on X)
    Xi = M'*(V{1}+D{1});
    for j = 2:(n_reg+1)
        Xi = Xi + V{j} + D{j};
    end
    X = IF*Xi;
    
    % solve the sparse noise term
    S = Y - V{1};
    S = sign(S).*max(abs(S) - lambda_2,0);

    % compute the Mourau proximity operators
    for j=1:(n_reg+1)
        %  data term (V1)
        if  reg(j) == 1
            V{j} = (1/(1+mu)*(Y-S+mu*(M*X-D{j})));
        end
        %  positivity   (V2)
        if  reg(j) == 2
            V{j} = max(X-D{j},0);
        end
        % addone  (project on the affine space sum(x) = 1)  (V3)
        if  reg(j) == 3
            nu_aux = X - D{j};
            V{j} = nu_aux + repmat((1-sum(nu_aux))/n,n,1);
        end
        % nuclear norm  (V4)
        if  reg(j) == 4
            V{j} = nuclear_norm_shrinkage(X-D{j},lambda_1/mu);
        end
        % TV  (V5 and V6)
        if  reg(j) == 5
            % update V5: solves the problem:
            %    min 0.5*||L*V5-(V6+D7)||^2+0.5*||V5-(X-d5)||^2
            %      V5
            %
            % update V6: min 0.5*||V6-(L*V5-D6)||^2 + lambda_tv * |||V6||_{1,1}
            
            nu_aux = X - D{j};
            % convert nu_aux into image planes
            % convert X into a cube
            nu_aux5_im = reshape(nu_aux',im_size(1), im_size(2), n);
            % compute V5 in the form of image planes
            for k =1:n
                % V5
                V5_im(:,:,k) = real(ifft2(IL.*fft2(DhH(V{j+1}{k}{1}+D{j+1}{k}{1}) ...
                    +  DvH(V{j+1}{k}{2}+D{j+1}{k}{2}) +  nu_aux5_im(:,:,k))));
                % V6
                aux_h = Dh(V5_im(:,:,k));
                aux_v = Dv(V5_im(:,:,k));
                if strcmp(tv_type, 'niso')  % non-isotropic TV
                    V{j+1}{k}{1} = soft(aux_h - D{j+1}{k}{1}, lambda_TV/mu);   %horizontal
                    V{j+1}{k}{2} = soft(aux_v - D{j+1}{k}{2}, lambda_TV/mu);   %vertical
                else    % isotropic TV
                    % Vectorial soft threshold
                    aux = max(sqrt((aux_h - D{j+1}{k}{1}).^2 + (aux_v - D{j+1}{k}{2}).^2)-lambda_TV/mu,0);
                    V{j+1}{k}{1} = aux./(aux+lambda_TV/mu).*(aux_h - D{j+1}{k}{1});
                    V{j+1}{k}{2} = aux./(aux+lambda_TV/mu).*(aux_v - D{j+1}{k}{2});
                end
                % update D6
                D{j+1}{k}{1} =  D{j+1}{k}{1} - (aux_h - V{j+1}{k}{1});
                D{j+1}{k}{2} =  D{j+1}{k}{2} - (aux_v - V{j+1}{k}{2});
            end
            % convert V6 to matrix format
            V{j} = reshape(V5_im, prod(im_size), n)';    
        end
    end
        
    
    % update Lagrange multipliers
    for j=1:(n_reg+1)
        if  reg(j) == 1
            D{j} = D{j} - (M*X-V{j});
        else
            D{j} = D{j} - (X-V{j});
        end
    end

    
    % compute residuals
    if mod(i,10) == 1
        st = [];
        for j = 1:(n_reg+1)
            if  reg(j) == 1
                res(j) = norm(M*X-V{j},'fro');
                st = strcat(st,sprintf(' res_p(%i) = %2.6f',reg(j),res(j) ));
            else
                res(j) = norm(X-V{j},'fro');
                st = strcat(st,sprintf(' res_p(%i) = %2.6f',reg(j),res(j) ));
            end
        end

        res_p = sum(abs(res));

        for j = 1:(n_reg+1)
            res(j) = mu*norm(V10{j}-V{j},'fro');
            st = strcat(st,sprintf(' res_d(%i) = %2.6f',reg(j),res(j) ));
        end

        res_d = sum(abs(res));

        if  strcmp(verbose,'yes')
            fprintf(strcat(sprintf('iter = %i -',i),st,'\n'));
        end
    
        % compute RMSE
        if true_x
            rmse(i)= norm(X-XT,'fro');
            if  strcmp(verbose,'yes')
                fprintf(strcat(sprintf('iter = %i - ||Xhat - X|| = %2.3f',i, rmse(i)),'\n'));
            end
            
        end
    end
    
    i = i+1;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
