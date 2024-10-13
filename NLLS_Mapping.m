% Description: T1, T2, PD and B0 mapping with nonlinear least-squares 
% estimation from phase-cycled bSSFP profiles. Also prodocues qualitative
% elliptical parameter maps

% This code is for research purposes only.

% Author of function: 
% Berk Can Acikgoz, Bern, Switzerland
% E-mail: berk.acikgoez@students.unibe.ch
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland


% profiles:    2D or 3D phase-cycled bSSFP profiles (Nx x Ny x Nz x NPc)
% TR:          Repetition time in miliseconds
% fa:          Excitation angle
% pc_step:     Step size between individual RF phase increments


%%% The outputs are quantitative T1,T2 and PD, and elliptical parameter 
%%% maps namely a, q and M

function [T1map, T2map, pdmap, b0map, amap, qmap, Mmap] = ...
    NLLS_Mapping(profiles, TR, fa, pc_step)


%%% Getting the correct dimensions for single-slice input
if length(size(profiles))==3
    profiles = permute(profiles, [1 2 4 3]);
end


%%% Initialize the maps with zeroes
T2map = zeros(size(profiles,1), size(profiles,2), size(profiles,3));
T1map = T2map;
amap = T2map;
qmap = T2map;
pdmap = T2map;
Mmap = T2map;
resmap = T2map;
b0map = T2map;
brain_mask = T2map;


%%% Intensity-based masking. Filters out voxels, based on z-score
%%% thresholding of magnitude of arithmetic mean of PC-bSSFP profiles
for k = 1:size(profiles,3)
    profiles_slice = squeeze(profiles(:,:,k,:));
    brain_mask(:,:,k) = IntensityMask(abs(sum(profiles_slice,3)), -0.5);
end



%%% NLLS fitting loop over voxels
for k = 1:size(profiles,3)

    profiles_slice = squeeze((profiles(:,:,k,:)));
    
    profiles_rot = profiles_slice.*repmat(exp(-1i*angle(mean(profiles_slice,3))), [1 1 size(profiles_slice,3)]);
    
    %%% Here, a dictionary of a single PC-bSSFP ellipse is created
    % In the dictionary, only the frequency offset df is varying. The 
    % dictionary is designed like this for phase-constrained least squares
    % fitting to make it well-conditioned. The outputs of the
    % PhaseCorrection function are phase corrected PC-bSSFP profiles and
    % the offresonance-induced phase map.  
    D = DictionaryCreation(TR, 80, 10, fa, pc_step);
    [phase_corrected_profiles, phases] = PhaseCorrection(profiles_slice, D);

    for i = 1:size(profiles_rot,1)
        for j = 1:size(profiles_slice,2)

            % NLLS fitting only if the voxel insensity is high enough to
            % gain time
            if brain_mask(i,j,k)>0

                

                % Converting phase of the voxel to off-resonance frequency
                % in Hz
                phase_hz = phases(i,j)/pi*(1e3/TR);

                % p is the single PC-bSSFP profile
                p = squeeze(phase_corrected_profiles(i,j,:));

                % Function call to NLLS-based parameter estimation based
                % on the elliptical model of the PC-bSSFP. The "sol" output
                % is the elliptical parameters a, q and M; quant is the
                % quantified parameters T1, T2, PD and B0; res is the
                % residual between the fitted and measured PC-bSSFP
                % profile.
                [sol, quant, res] = ...
                    bSSFP_NonlinLS_Fit_WithPhase_Quant(...
                                            p, phase_hz, pc_step, TR, fa);
            
                T1map(i,j,k) = quant(1);
                T2map(i,j,k) = quant(2);
                pdmap(i,j,k) = quant(3);  

                amap(i,j,k) = sol(1);
                qmap(i,j,k) = sol(2);
                Mmap(i,j,k) = sol(3);
               
                resmap(i,j,k) = res;
            
            else
                T2map(i,j,k) = 0;
                pdmap(i,j,k) = 0;
                T1map(i,j,k) = 0;
                resmap(i,j,k) = 0;
                qmap(i,j,k) = 0;
                Mmap(i,j,k) = 0;
                amap(i,j,k) = 0;
            end
        end
        b0map(:,:,k) = phases;

    end
    progressPercent = (k/size(profiles,3))*100;
    fprintf('NLLS Fitting Progress: %.2f%%\n', progressPercent);
    end
    
    T2map = squeeze(abs(T2map));
    T1map = squeeze(abs(T1map));
    amap = squeeze(abs(amap));
    qmap = squeeze(abs(qmap));
    Mmap = squeeze(abs(Mmap));
    pdmap = squeeze(abs(pdmap));

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
    
function e1 = getE1(a,q, alpha)
    %%% This function calculates E1 based on the elliptical parameters a
    %%% and q, and flip angle alpha
    c = cosd(alpha);
    e1 = (-q-c.*(q.*(a.^2))+a+(a.*c))./(a+(a.*c)-(q.*c)-(a.^2).*q);
end

function [intense_mask, intense_mask_nand] = IntensityMask(profiles,threshold,op)

    %%% This function masks the slices based on the magnitude of the
    %%% complex sum of PC-bSSFP profiles. This is used only to speed the
    %%% NLLS process up by excluding background voxels.
    complex_sum = abs(sum(profiles,3));
    
    %%& Mean and standard deviation of the magnitude of the complex sum
    % are set to 0 and 1, respectively, by standardizing the image.
    complex_sum_feature = StandardizeImage(complex_sum);
    
    %%% Regions below a certain threshold is set to bo 0. Usually -0.7 as a
    %%% threshold works very well.
    intense_mask = 1*(complex_sum_feature>threshold);

    %%% Another version of the mask that has NaNs instead of zeros, might
    %%% come in handy for unwrapping etc.
    intense_mask_nand = 1*intense_mask;
    intense_mask_nand(intense_mask==0) = nan;

end




function data = StandardizeImage(data)
    %%% Accepts an image, normalizes it such that its mean is 0 and its
    %%% standard deviation is 1.
    sizerarr = size(data);
    data = data(:);
    data = data-mean(data);
    data = data/(std(data));
    data = reshape(data, sizerarr);
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


function [sol, quant,res,found_profile] = bSSFP_NonlinLS_Fit_WithPhase_Quant(p, phase, pc_step, tr, fa)
    
    if nargin<3
        pc_step = 360/length(p);
    end
    
    pc = (0:pc_step:359)';
    pc = pc*pi/180;
    phase = phase*2*pi*(tr/1e3);
    pr = phase/2;
    
    %%% This inline function is an implementation of canonical elliptical
    %%% equation where the off-resonance is also fixed ("phase" in the 
    % function)
    ellipse_canon = @(x) (x(3)*(((1-x(1)*exp(1i*(-phase+pc))))...
        ./(1-x(2)*cos((-phase+pc)))))*exp(1i*pr);

    %%% Difference function between the fitted and measured PC-bSSFP
    %%% profiles
    f = @(x) reimconcat( -ellipse_canon(x) ) + reimconcat(p);

    %%% Initial guess as a=0.9, q=0.1, M=1
    u0 = [ 0.9; 0.1; 1];

    %%% This is the term that accounts for the rotation of the ellipse
    % around the origin (not around itself), note that the initial
    %  "ellipse_cannon" implements only the PC-bSSFP rotation around
    % itself, not around the origin. It is fixed with the phase-kern
    phase_kern = [cos(pr)*eye(length(pc)) -sin(pr)*eye(length(pc));...
        sin(pr)*eye(length(pc)) cos(pr)*eye(length(pc))];

    %%% Partial derivative calculation functions of the difference function
    % with respect to the 3 elliptical parameters 
    j1c = @(x) x(3)*(exp(1i*(-phase+pc))./(1-x(2)*cos(-phase+pc)));
    J1 = @(x) phase_kern*[real(j1c(x)); imag(j1c(x))];

    j2c = @(x) -x(3)* ( (1-x(1)*exp(1i*(-phase+pc))).*cos(-phase+pc) ./ ( ( 1-x(2)*cos(-phase+pc) ).^2 ) );
    J2 = @(x) phase_kern*[real(j2c(x)); imag(j2c(x))];


    j3c = @(x) -1*(((1-x(1)*exp(1i*(-phase+pc))))./(1-x(2)*cos(-phase+pc))); 
    J3 = @(x) phase_kern*[real(j3c(x)); imag(j3c(x))];

    %%% Jacobian and cost function calculation for the iterations of
    %%% Guass-Newton
    J = @(x) [J1(x)'; J2(x)'; J3(x)']';
    jf_inline = @(x) jf(x,J,f);

    u = u0;

    for i = 1:50

        [Jm, F] = jf_inline(u);
    
        %%% Gauss-Newton update of the current solution
        h = pinv(Jm)*(-F);
        un = u + h;
        
        %%% If the change is insignificant, stop the iterations 
        % (local minimum)
        if (norm(un-u,2)/norm(u,2))<1e-3
            break
        end

        u = un(1:3,1);
    end
        
        %%% Preparing the outputs based on the fitted elliptical parameters
        sol = u;
        sol = abs(sol);
        found_profile = ellipse_canon(sol);
        res = norm(f(sol),2);
        a = abs(sol(1));
        q = abs(sol(2));
        M = abs(sol(3));
        T2 = tr/log(1/a);
        E1 = abs(getE1(a,q,fa));
        T1 = abs(tr/log(1/E1));
        pd = (M/((1-E1)*sind(fa)/(1-E1*cosd(fa)-(sol(1)^2)*(E1-cosd(fa)))))...
            /(exp(-(tr/2)/T2));
        quant = [T1, T2, pd];
end

function [Jm, F] = jf(x,J,f)

    Jm = J(x);
    F = f(x);

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

function concat_res = reimconcat(signal, dim)

    if nargin == 1
        signal_re = real(signal);
        signal_im = imag(signal);
        concat_res = [signal_re;signal_im];
    else
        signal_re = real(signal);
        signal_im = imag(signal);
        concat_res = cat(dim, signal_re, signal_im);
    end

end

