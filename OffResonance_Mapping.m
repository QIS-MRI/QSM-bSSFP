% Description: Off resonance mapping for phase-cycled bSSFP

% This code is for research purposes only.

% Author of function: 
% Berk Can Acikgoz, Bern, Switzerland
% E-mail: berk.acikgoez@students.unibe.ch
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland


% profiles:    2D or 3D phase-cycled bSSFP profiles (Nx x Ny x Nz x NPc)
% TR:          Repetition time in miliseconds

function [offres] = OffResonance_Mapping(profiles, TR)
    %%% Getting the correct dimensions for single-slice input
    if length(size(profiles))==3
        profiles = permute(profiles, [1 2 4 3]);
    end

    %%% Rotating profiles onto the real axis
    profiles = profiles.*repmat(exp(-1i*angle(mean(profiles,4))), [1 1 1 18]);

    %%% Here, a dictionary of a single PC-bSSFP ellipse is created
    % In the dictionary, only the frequency offset df is varying. The 
    % dictionary is designed like this for phase-constrained least squares
    % fitting to make it well-conditioned. The outputs of the
    % PhaseCorrection function are phase corrected PC-bSSFP profiles and
    % the offresonance-induced phase map.  
    D = DictionaryCreation(TR, 80, 10, 15, 360/size(profiles,4));

    for slice = 1:size(profiles, 3)
        slice_profiles = squeeze(profiles(:,:,slice,:));
        [phase_corrected_profiles(:,:,slice,:), phases(:,:,slice)] = ...
            PhaseCorrection(slice_profiles, D);
        progressPercent = 100*slice/size(profiles,3);
        fprintf('Offresonance Mapping Progress: %.2f%%\n', progressPercent);
        
    end

    phases_dict = phases;
    offres = phases_dict;  

end
function [D, parameters] = DictionaryCreation(tr,t2arr,rtrarr,fa, pc_step)
    
    %%% This function creates a dictionary based on given T2 and T1/T2 (RTR)
    %%% ranges. Steps in off-resonance is equal to the number of
    %%% phase-cycles. It uses the analytic signal equation of PC-bSSFP
    te = tr/2;
    freq_step = 1e3/tr/(360/pc_step);
    [b0, rtr, t2] = ndgrid(linspace(0,1e3/tr-freq_step,(360/pc_step)),rtrarr, t2arr);
    D = 1i+zeros(round(360/pc_step), numel(rtr));
    
    for i = 1:numel(rtr)
        p = bSSFPAnalytic(rtr(i)*t2(i), t2(i), te, tr, fa, 0:pc_step:359, b0(i));
        D(:,i) = ((squeeze((p))));
    end
    parameters = [b0(:), rtr(:), t2(:)];

end

function [corrected_profiles, phase_offset] = PhaseCorrection(profiles, dictionary)

    %%% This function implements the phase correction algorithm presented
    % in the paper. It is a real-constrained least-squares problem where
    % solution weights are also constrained to share a constant phase. Its
    % solution is inspired by:
    % Bydder, M. (2010). Solution of a complex least squares problem with 
    % constrained phase. Linear algebra and its applications, 
    % 433(11-12), 1719-1721.

    D = dictionary;

    %%% Precalculate the M for speed
    M = pinv(real(D'*D));

    %%% Initialize the outputs
    corrected_profiles =1i+zeros(size(profiles));
    phase_offset = zeros(size(profiles,1), size(profiles,2));

    for x = 1:size(profiles,1)
        for y = 1:size(profiles,2)

            p = squeeze(profiles(x,y,:));
            h = D'*p;
            c = h.'*M*h;

            b = angle(c);
            a = b/2;

            cp = p*exp(-1i*a);
            corrected_profiles(x,y,:) = cp;

            phase_offset(x,y) = -a;

        end
    end

end

function signal = bSSFPAnalytic(t1,t2,te,tr,alpha,pc,off_res)
    %%% PC-bSSFP profile simulation with analytic equation
    e1 = exp(-tr/t1);
    e2 = exp(-tr/t2);
    
    alpha = alpha*pi/180;
    d = (1-e1*cos(alpha)-(e2^2)*(e1-cos(alpha)));
    b = e2*(1-e1)*(1+cos(alpha))/d;
    a = e2;
    M = (1-e1)*sin(alpha)/d;
    delt = pc*pi/180 ;
    theta0 = 2*pi*off_res*(tr*1e-3);
    theta = theta0-delt;
    phi = theta0*te/tr;
    signal = (M*(1-a*exp(1i*theta))./(1-b*cos(theta)))*exp(-1i*phi);
    signal = conj(signal)*exp(-te/t2);

end

function [a] = FitEllipseStable(x, y)
    
    %%% Stable elliptical fitting, implementation of 
    % Halır, R., & Flusser, J. (1998, February). Numerically stable direct
    % least squares fitting of ellipses. In Proc. 6th International 
    % Conference in Central Europe on Computer Graphics and Visualization.
    % WSCG (Vol. 98, pp. 125-132). Plzen-Bory: Citeseer.
    % Directly implemented as explained in the original work.
    
    D1 = [x .^ 2, x .* y, y .^ 2]; % quadratic part of the design matrix
    D2 = [x, y, ones(size(x))]; % linear part of the design matrix
    D1 = D1; D2 = D2;
    S1 = D1' * D1; % quadratic part of the scatter matrix
    S2 = D1' * D2; % combined part of the scatter matrix
    S3 = D2' * D2; % linear part of the scatter matrix
    T = - inv(S3) * S2'; % for getting a2 from a1
    M = S1 + S2 * T; % reduced scatter matrix   
    M = [M(3, :) ./ 2; - M(2, :); M(1, :) ./ 2]; % premultiply by inv(C1)
    [evec, eval] = eig(M); % solve eigensystem
    cond = 4 * evec(1, :) .* evec(3, :) - evec(2, :) .^ 2; % evaluate a^Ca
    c = find(cond > 0);
    if isempty(c)
        a1 = ones(3,1);
    else
        a1 = evec(:, c(1)); % eigenvector for min. pos. eigenvalue
    end
    a = [a1; T * a1];

end