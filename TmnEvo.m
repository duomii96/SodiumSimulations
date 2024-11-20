% % % % % % % % % % % % % %
% % Tmn class 
% %  Simulation of Na Dynamics
% %  Hard Pulses and Relaxation with or without a macroscopic anisotropy
% %  and a potential B0shift/inhomogeneity
% %  
% %  Tmn matrix contains the the fraction of each Tmn basis element

classdef TmnEvo < handle
properties(Constant)
    Teq = [0 0 0 1 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0]; % thermal equilibrium = T10
end
properties
    
    Tmn = [0 0 0 1 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0];
    % % % % % % % % % % % % % % % % % % % % % 
    % % order m = 1, 2, 3;                   --> changes for relaxation, constant for pulses
    % % rank  n = -3, -2, -1, 0, 1 , 2, 3    --> changes for pulses, constant for relaxation
    % % Tmn  = [T1-3 T1-2 ... T13;
    % %         T2-3 T2-2 ... T23; 
    % %         T3-3 T3-2 ... T33; ]
    % % % % % % % % % % % % % % % % % % % % %
    
    B0             = 9.4;       % T
    w0                          % Lamor frequency
    tauC                        % correlation time
    wQ                          % wQ_RMS 
    wQbar          = 0;         % mean(wQ(t)) --> anisotropy
    Jen                         % Jen value for Jen model
    tsDataPoints   = 2048;      % Number of points for evolution period (old 2048)
    
    
    wShift         = 0;         % B-field inhomogeneity -> mean(Omega)
    wShift_RMS     = 0;         % B-field inhomogeneity -> ~std(Omega)  --> intrinsical, no refocus by 180
    
    wShift_FID     = 0;         % B-field inhomogeneity for 1. Dim FID  --> refocus by 180 possible (but not performed in FID
    wShift_Freq    = 0;  
    
    DeltaHz        = 0;         % off-resonance frequency, Hz
    noFrames       = 1e3;       % number of Frames for Soft Pulse calculation


end
methods(Static)  
    % returns Tmn amplitude for given Tmn Matrix
    function [Tmn_value] = getTmn_matrix(Tmn, m, n)
        Tmn_value = Tmn(m, n+4);
    end
    
    % convert symmetric/asymmetric Tmn to Tmn with plus/minus
    function [Tmn_a, Tmn_s] = TmnPM2TmnAS(Tmn_m, Tmn_p)
        Tmn_a = 1/sqrt(2)*(Tmn_m-Tmn_p);
        Tmn_s = 1/sqrt(2)*(Tmn_m+Tmn_p);
    end
    
    % convert Tmn with plus/minus to symmetric/asymmetric Tmn
    function [Tmn_m, Tmn_p] = TmnAS2TmnPM(Tmn_a, Tmn_s)
        Tmn_p = 1/sqrt(2)*(Tmn_s-Tmn_a);
        Tmn_m = 1/sqrt(2)*(Tmn_s+Tmn_a);
    end
    
    % converts M vector to Tmn matrix
    function [Tmn] = M2Tmn(M)
      	T10         = M(1);
        [T1m1, T11] = TmnEvo.TmnAS2TmnPM(M(2), M(3));
        
        T20         = M(4);
        [T2m1, T21] = TmnEvo.TmnAS2TmnPM(M(5), M(6));
        [T2m2, T22] = TmnEvo.TmnAS2TmnPM(M(7), M(8));
        
        T30         = M(9);
        [T3m1, T31] = TmnEvo.TmnAS2TmnPM(M(10), M(11));
        [T3m2, T32] = TmnEvo.TmnAS2TmnPM(M(12), M(13));
        [T3m3, T33] = TmnEvo.TmnAS2TmnPM(M(14), M(15));
        
        Tmn = [0 0 T1m1 T10 T11 0 0; 0 T2m2 T2m1 T20 T21 T22 0; T3m3 T3m2 T3m1 T30 T31 T32 T33];
    end
    
    % converts M vector timeseries of form M(time dim, Tmn dim) to Tmn matrix time series
    function [Tmns] = Ms2Tmns(M)
        sizeM       = size(M); 
        Tmn_non     = zeros([1 sizeM(1)]); % nonexisting Tmns
        
      	T10         = M(:,1);
        [T1m1, T11] = TmnEvo.TmnAS2TmnPM(M(:,2), M(:,3));
        
        T20         = M(:,4);
        [T2m1, T21] = TmnEvo.TmnAS2TmnPM(M(:,5), M(:,6));
        [T2m2, T22] = TmnEvo.TmnAS2TmnPM(M(:,7), M(:,8));
        
        T30         = M(:,9);
        [T3m1, T31] = TmnEvo.TmnAS2TmnPM(M(:,10), M(:,11));
        [T3m2, T32] = TmnEvo.TmnAS2TmnPM(M(:,12), M(:,13));
        [T3m3, T33] = TmnEvo.TmnAS2TmnPM(M(:,14), M(:,15));
        
        Tmns = zeros([3 7 sizeM(1)]);
        Tmns(1,:,:) = [ Tmn_non; Tmn_non; T1m1'; T10'; T11'; Tmn_non; Tmn_non];
        Tmns(2,:,:) = [ Tmn_non;   T2m2'; T2m1'; T20'; T21';    T22'; Tmn_non];
        Tmns(3,:,:) = [   T3m3';   T3m2'; T3m1'; T30'; T31';    T32';    T33'];
        
    end   
    
    
    % converts Tmn matrix to M vector
    function [M] = Tmn2M(Tmn)
         M             = zeros(15,1);
      	 M(1)          = TmnEvo.getTmn_matrix(Tmn, 1, 0);
        [M(2), M(3)]   = TmnEvo.TmnPM2TmnAS(TmnEvo.getTmn_matrix(Tmn, 1, -1), TmnEvo.getTmn_matrix(Tmn, 1, 1));
        
         M(4)          = TmnEvo.getTmn_matrix(Tmn, 2, 0);
        [M(5), M(6)]   = TmnEvo.TmnPM2TmnAS(TmnEvo.getTmn_matrix(Tmn, 2, -1), TmnEvo.getTmn_matrix(Tmn, 2, 1));
        [M(7), M(8)]   = TmnEvo.TmnPM2TmnAS(TmnEvo.getTmn_matrix(Tmn, 2, -2), TmnEvo.getTmn_matrix(Tmn, 2, 2));
        
         M(9)          = TmnEvo.getTmn_matrix(Tmn, 3, 0);
        [M(10), M(11)] = TmnEvo.TmnPM2TmnAS(TmnEvo.getTmn_matrix(Tmn, 3, -1), TmnEvo.getTmn_matrix(Tmn, 3, 1));
        [M(12), M(13)] = TmnEvo.TmnPM2TmnAS(TmnEvo.getTmn_matrix(Tmn, 3, -2), TmnEvo.getTmn_matrix(Tmn, 3, 2));
        [M(14), M(15)] = TmnEvo.TmnPM2TmnAS(TmnEvo.getTmn_matrix(Tmn, 3, -3), TmnEvo.getTmn_matrix(Tmn, 3, 3));
    end

end
methods
    % constructor
    function obj = TmnEvo(B0, tauC, wQ, wQbar, Jen, wShift, wShift_RMS, wShift_FID, Tmn_value)%, wShift_Freq
        if nargin == 8 % default is eq matrix = T10
            obj.Tmn = TmnEvo.Teq;
        else
            obj.Tmn = Tmn_value;
        end
        obj.B0         = B0;
        obj.w0         = getValues.get_w0(B0);
        
        obj.tauC       = tauC;
        obj.wQ         = wQ;
        obj.wQbar      = wQbar;
        obj.Jen        = Jen;
        
        obj.wShift     = wShift;
        obj.wShift_RMS = wShift_RMS;
        obj.wShift_FID = wShift_FID;
%         obj.wShift_Freq = wShift_Freq;
    end
    
    % % copy TmnEvo oject
    function [Tmn] = copy(obj)
        Tmn = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift, obj.wShift_RMS, obj.wShift_FID, obj.Tmn);
        Tmn.tsDataPoints = obj.tsDataPoints; 
    end
    
    % % getter and setter functions for a single Tmn element
    function [Tmn_value] = getTmn(obj, m, n)
        Tmn_value = obj.Tmn(m, n+4);
    end
    function setTmn(obj, m, n, Tmn_value)
        obj.Tmn(m, n+4) = Tmn_value;
    end
    

    
    % % return TmnEvo object after hard pulse application
    function [Tmn_pulse] = pulse(obj, flipAngle, phase)
        Tmn_prePulse = obj.Tmn;
        % % for all m Tmn,new = (P * Tmn^T)^T      ^T = transpose operator
%         getValues.get_pulseOperator(flipAngle, phase, 1)
        Tmn_pulse1 = (getValues.get_pulseOperator(flipAngle, phase, 1) * Tmn_prePulse(1,:).').';
        Tmn_pulse2 = (getValues.get_pulseOperator(flipAngle, phase, 2) * Tmn_prePulse(2,:).').';
        Tmn_pulse3 = (getValues.get_pulseOperator(flipAngle, phase, 3) * Tmn_prePulse(3,:).').';
        Tmn_pulse_matrix  = [Tmn_pulse1 ; Tmn_pulse2; Tmn_pulse3]  ; 
        Tmn_pulse = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift, obj.wShift_RMS, obj.wShift_FID, Tmn_pulse_matrix);
    end
    
%     % % relaxation 
%     % % without anisotropy, B0 shift, etc
%     % % return: Tmn object after relaxation, 
%     % %         Tmn Matrix with time evolution
%     function [Tmn_relaxation, Tmn_timeSeries] = relaxation(obj, tevo)
%         Tmn_preRelax = obj.copy(); 
%         Tmn_preRelax.Tmn = Tmn_preRelax.Tmn - Tmn_preRelax.Teq; % recover thermal equilibrium
%         
%         ts = 0:tevo/obj.tsDataPoints:tevo;
%         Tmn_non = zeros([1 length(ts)]); % nonexisting Tmns
%         T10_ts = ones([1 length(ts)]);   
%                 
%         fm333  = getValues.get_f(-3, ts, obj.tauC, obj.wQ, obj.w0);
%         fm2kks = getValues.get_f(-2, ts, obj.tauC, obj.wQ, obj.w0);
%         fm1kks = getValues.get_f(-1, ts, obj.tauC, obj.wQ, obj.w0);
%         f0kks  = getValues.get_f( 0, ts, obj.tauC, obj.wQ, obj.w0);
%         f1kks  = getValues.get_f( 1, ts, obj.tauC, obj.wQ, obj.w0);
%         f2kks  = getValues.get_f( 2, ts, obj.tauC, obj.wQ, obj.w0);        
%         f333   = getValues.get_f( 3, ts, obj.tauC, obj.wQ, obj.w0);
%           
%         % zero quantum relaxation
%         f011 = f0kks(1,:); f013 = f0kks(2,:); f033 = f0kks(3,:); f022 = f0kks(4,:);
%         T10_pre = Tmn_preRelax.getTmn(1, 0);
%         T30_pre = Tmn_preRelax.getTmn(3, 0);
%         T20_pre = Tmn_preRelax.getTmn(2, 0);
%         T10_relax = f011*T10_pre+f013*T30_pre + T10_ts;  % add T10_ts to recover thermal equilibrium
%         T30_relax = f013*T10_pre+f033*T30_pre;
%         T20_relax = f022* T20_pre;
%         
%         % single quantum relaxation
%         fm111 = fm1kks(1,:); fm113 = fm1kks(2,:); fm133 = fm1kks(3,:); fm122 = fm1kks(4,:);
%         f111 = f1kks(1,:); f113 = f1kks(2,:); f133 = f1kks(3,:); f122 = f1kks(4,:);
%         T1m1_pre = Tmn_preRelax.getTmn(1, -1); T11_pre = Tmn_preRelax.getTmn(1, 1);
%         T3m1_pre = Tmn_preRelax.getTmn(3, -1); T31_pre = Tmn_preRelax.getTmn(3, 1);
%         T2m1_pre = Tmn_preRelax.getTmn(2, -1); T21_pre = Tmn_preRelax.getTmn(2, 1);
%         % n=-1
%         T1m1_relax = fm111*T1m1_pre+fm113*T3m1_pre;   
%         T3m1_relax = fm113*T1m1_pre+fm133*T3m1_pre;
%         T2m1_relax = fm122* T2m1_pre;
%         % n=+1
%         T11_relax = f111*T11_pre+f113*T31_pre;   
%         T31_relax = f113*T11_pre+f133*T31_pre;
%         T21_relax = f122* T21_pre;
%         
%         % multiple quantum relaxation
%         % m = 2
%         f222 = f2kks(1,:); f233 = f2kks(2,:);
%         fm222 = fm2kks(1,:); fm233 = fm2kks(2,:);
%         T2m2_pre = Tmn_preRelax.getTmn(2, -2); T22_pre = Tmn_preRelax.getTmn(2, 2);
%         T3m2_pre = Tmn_preRelax.getTmn(3, -2); T32_pre = Tmn_preRelax.getTmn(3, 2);
%         % n=-2
%         T2m2_relax = fm222* T2m2_pre;
%         T3m2_relax = fm233* T3m2_pre;
%         % n=+2
%         T22_relax = f222* T22_pre;
%         T32_relax = f233* T32_pre; 
%         % m = 3
%         T3m3_pre = Tmn_preRelax.getTmn(3, -3); T33_pre = Tmn_preRelax.getTmn(3, 3);
%         T3m3_relax = fm333* T3m3_pre;
%         T33_relax  = f333* T33_pre;
%         
%         Tmn_timeSeries = zeros([3 7 length(ts)]);
%         Tmn_timeSeries(1,:,:) = [   Tmn_non;    Tmn_non; T1m1_relax; T10_relax; T11_relax;   Tmn_non;   Tmn_non];
%         Tmn_timeSeries(2,:,:) = [   Tmn_non; T2m2_relax; T2m1_relax; T20_relax; T21_relax; T22_relax;   Tmn_non];
%         Tmn_timeSeries(3,:,:) = [T3m3_relax; T3m2_relax; T3m1_relax; T30_relax; T31_relax; T32_relax; T33_relax];
%         
%         Tmn_relax = Tmn_timeSeries(:,:,end);        
%         Tmn_relaxation = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift, obj.wShift_RMS, Tmn_relax);
%     end  
       
%     % % relaxation  with macroscopic anisotropy, without B0 inhomogeneity
%     % % return: Tmn object after relaxation, 
%     % %         Tmn Matrix with time evolution
%     function [Tmn_relaxation, Tmn_timeSeries] = relaxation_aniso(obj, tevo)
%         Tmn_preRelax = obj.copy(); 
%         Tmn_preRelax.Tmn = Tmn_preRelax.Tmn - Tmn_preRelax.Teq;  % recover thermal equilibrium
% 
%         ts = 0:tevo/obj.tsDataPoints:tevo;
%         Tmn_non = zeros([1 length(ts)]); % nonexisting Tmns
%         T10_ts = ones([1 length(ts)]);   
%         
%         fm333  = getValues.get_f_aniso(-3, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         fm2kks = getValues.get_f_aniso(-2, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         fm1kks = getValues.get_f_aniso(-1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         f0kks  = getValues.get_f_aniso( 0, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         f1kks  = getValues.get_f_aniso( 1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         f2kks  = getValues.get_f_aniso( 2, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         f333   = getValues.get_f_aniso( 3, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%           
%         % zero quantum relaxation 
%         f011 = f0kks(1,:); f013 = f0kks(2,:); f033 = f0kks(3,:); f022 = f0kks(4,:);
%         T10_pre = Tmn_preRelax.getTmn(1, 0);
%         T30_pre = Tmn_preRelax.getTmn(3, 0);
%         T20_pre = Tmn_preRelax.getTmn(2, 0);
%         T10_relax = f011*T10_pre+f013*T30_pre + T10_ts;   % add T10_ts to recover thermal equilibrium
%         T30_relax = f013*T10_pre+f033*T30_pre;
%         T20_relax = f022* T20_pre;
%         
%         % single quantum relaxation
%         fm111 = fm1kks(1,:); fm113 = fm1kks(2,:); fm133 = fm1kks(3,:); fm122 = fm1kks(6,:);
%         fm112 = -fm1kks(4,:); fm123 = -fm1kks(5,:);
%         f111 = f1kks(1,:); f113 = f1kks(2,:); f133 = f1kks(3,:); f122 = f1kks(4,:);
%         f112 = f1kks(4,:); f123 = f1kks(5,:);
%         T1m1_pre = Tmn_preRelax.getTmn(1, -1); T11_pre = Tmn_preRelax.getTmn(1, 1);
%         T3m1_pre = Tmn_preRelax.getTmn(3, -1); T31_pre = Tmn_preRelax.getTmn(3, 1);
%         T2m1_pre = Tmn_preRelax.getTmn(2, -1); T21_pre = Tmn_preRelax.getTmn(2, 1);
%         % n=-1
%         T1m1_relax = fm111*T1m1_pre + fm112*T2m1_pre + fm113*T3m1_pre;   
%         T3m1_relax = fm113*T1m1_pre + fm123*T2m1_pre + fm133*T3m1_pre;
%         T2m1_relax = fm112*T1m1_pre + fm122*T2m1_pre + fm123*T3m1_pre;
%         % n=+1
%         T11_relax = f111*T11_pre + f112*T21_pre + f113*T31_pre;   
%         T31_relax = f113*T11_pre + f123*T21_pre + f133*T31_pre;
%         T21_relax = f112*T11_pre + f122*T21_pre + f123*T31_pre;
%         
%         % multiple quantum relaxation
%         % m = 2
%         f222 = f2kks(1,:); f233 = f2kks(2,:); f223 = f2kks(3,:); 
%         fm222 = fm2kks(1,:); fm233 = fm2kks(2,:); fm223 = -f2kks(3,:);
%         T2m2_pre = Tmn_preRelax.getTmn(2, -2); T22_pre = Tmn_preRelax.getTmn(2, 2);
%         T3m2_pre = Tmn_preRelax.getTmn(3, -2); T32_pre = Tmn_preRelax.getTmn(3, 2);
%         % n=-2
%         T2m2_relax = fm222*T2m2_pre + fm223*T3m2_pre;
%         T3m2_relax = fm233*T3m2_pre + fm223*T2m2_pre;
%         % n=+2
%         T22_relax = f222*T22_pre + f223*T32_pre;
%         T32_relax = f233*T32_pre + f223*T22_pre;
%         % m = 3 
%         T3m3_pre = Tmn_preRelax.getTmn(3, -3); T33_pre = Tmn_preRelax.getTmn(3, 3);
%         T3m3_relax = fm333* T3m3_pre;
%         T33_relax  = f333* T33_pre;
%         
%         Tmn_timeSeries = zeros([3 7 length(ts)]);
%         Tmn_timeSeries(1,:,:) = [   Tmn_non;    Tmn_non; T1m1_relax; T10_relax; T11_relax;   Tmn_non;   Tmn_non];
%         Tmn_timeSeries(2,:,:) = [   Tmn_non; T2m2_relax; T2m1_relax; T20_relax; T21_relax; T22_relax;   Tmn_non];
%         Tmn_timeSeries(3,:,:) = [T3m3_relax; T3m2_relax; T3m1_relax; T30_relax; T31_relax; T32_relax; T33_relax];
%         
%         Tmn_relax = Tmn_timeSeries(:,:,end);        
%         Tmn_relaxation = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift, obj.wShift_RMS, Tmn_relax);
%     end 
 


%     % % relaxation  with macroscopic anisotropy and chemical shift (Cauchy/Lorentzian distribution for offset)
%     % % return: Tmn object after relaxation, 
%     % %         Tmn Matrix with time evolution    
%      function [Tmn_relaxation, Tmn_timeSeries] = relaxation_shift(obj, tevo)
%         Tmn_preRelax = obj.copy(); 
%         Tmn_preRelax.Tmn = Tmn_preRelax.Tmn - Tmn_preRelax.Teq; % recover thermal equilibrium
%  
%         ts = 0:tevo/obj.tsDataPoints:tevo;
%         Tmn_non = zeros([1 length(ts)]); % nonexisting Tmns
%         T10_ts = ones([1 length(ts)]);  
%   
%         fm333  = getValues.get_f_aniso(-3, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         fm2kks = getValues.get_f_aniso(-2, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         fm1kks = getValues.get_f_aniso(-1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         f0kks  = getValues.get_f_aniso( 0, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         f1kks  = getValues.get_f_aniso( 1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         f2kks  = getValues.get_f_aniso( 2, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%         f333   = getValues.get_f_aniso( 3, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
%           
%         % zero quantum relaxation 
%         f011 = f0kks(1,:); f013 = f0kks(2,:); f033 = f0kks(3,:); f022 = f0kks(4,:);
%         T10_pre = Tmn_preRelax.getTmn(1, 0);
%         T30_pre = Tmn_preRelax.getTmn(3, 0);
%         T20_pre = Tmn_preRelax.getTmn(2, 0);
%         % T10_30_relax = [f011 f013; f013; f033]*[T10_pre; T30_pre];
%         % T10_relax = T10_30_relax(1);   T30_relax = T10_30_relax(2);
%         T10_relax = f011*T10_pre+f013*T30_pre + T10_ts;   % add T10_ts to recover thermal equilibrium
%         T30_relax = f013*T10_pre+f033*T30_pre;
%         T20_relax = f022* T20_pre;
%         
%         % single quantum relaxation
%         fm111 = fm1kks(1,:); fm113 = fm1kks(2,:); fm133 = fm1kks(3,:); fm122 = fm1kks(6,:);
%         fm112 = -fm1kks(4,:); fm123 = -fm1kks(5,:);
%         f111 = f1kks(1,:); f113 = f1kks(2,:); f133 = f1kks(3,:); f122 = f1kks(4,:);
%         f112 = f1kks(4,:); f123 = f1kks(5,:);
%         T1m1_pre = Tmn_preRelax.getTmn(1, -1); T11_pre = Tmn_preRelax.getTmn(1, 1);
%         T3m1_pre = Tmn_preRelax.getTmn(3, -1); T31_pre = Tmn_preRelax.getTmn(3, 1);
%         T2m1_pre = Tmn_preRelax.getTmn(2, -1); T21_pre = Tmn_preRelax.getTmn(2, 1);
%         % n=-1
%         T1m1_relax = (fm111*T1m1_pre + fm112*T2m1_pre + fm113*T3m1_pre).*exp(-1*(obj.wShift_RMS+1i*obj.wShift)*ts);   
%         T3m1_relax = (fm113*T1m1_pre + fm123*T2m1_pre + fm133*T3m1_pre).*exp(-1*(obj.wShift_RMS+1i*obj.wShift)*ts);
%         T2m1_relax = (fm112*T1m1_pre + fm122*T2m1_pre + fm123*T3m1_pre).*exp(-1*(obj.wShift_RMS+1i*obj.wShift)*ts);
%         % n=+1
%         T11_relax = (f111*T11_pre + f112*T21_pre + f113*T31_pre).*exp(-1*(obj.wShift_RMS-1i*obj.wShift)*ts);   
%         T31_relax = (f113*T11_pre + f123*T21_pre + f133*T31_pre).*exp(-1*(obj.wShift_RMS-1i*obj.wShift)*ts);
%         T21_relax = (f112*T11_pre + f122*T21_pre + f123*T31_pre).*exp(-1*(obj.wShift_RMS-1i*obj.wShift)*ts);
%         
%         % multiple quantum relaxation
%         % m = 2
%         f222 = f2kks(1,:); f233 = f2kks(2,:); f223 = f2kks(3,:); 
%         fm222 = fm2kks(1,:); fm233 = fm2kks(2,:); fm223 = f2kks(3,:);
%         T2m2_pre = Tmn_preRelax.getTmn(2, -2); T22_pre = Tmn_preRelax.getTmn(2, 2);
%         T3m2_pre = Tmn_preRelax.getTmn(3, -2); T32_pre = Tmn_preRelax.getTmn(3, 2);
%         % n=-2
%         T2m2_relax = (fm222*T2m2_pre + fm223*T3m2_pre).*exp(-2*(obj.wShift_RMS+1i*obj.wShift)*ts);
%         T3m2_relax = (fm233*T3m2_pre + fm223*T2m2_pre).*exp(-2*(obj.wShift_RMS+1i*obj.wShift)*ts);
%         % n=+2
%         T22_relax = (f222*T22_pre + f223*T32_pre).*exp(-2*(obj.wShift_RMS-1i*obj.wShift)*ts);
%         T32_relax = (f233*T32_pre + f223*T22_pre).*exp(-2*(obj.wShift_RMS-1i*obj.wShift)*ts);
%         % m = 3 
%         T3m3_pre = Tmn_preRelax.getTmn(3, -3); T33_pre = Tmn_preRelax.getTmn(3, 3);
%         T3m3_relax = fm333* T3m3_pre.*exp(-3*(obj.wShift_RMS+1i*obj.wShift)*ts);
%         T33_relax  = f333* T33_pre.*exp(-3*(obj.wShift_RMS-1i*obj.wShift)*ts);
% 
%         Tmn_timeSeries = zeros([3 7 length(ts)]);
%         Tmn_timeSeries(1,:,:) = [   Tmn_non;    Tmn_non; T1m1_relax; T10_relax; T11_relax;   Tmn_non;   Tmn_non];
%         Tmn_timeSeries(2,:,:) = [   Tmn_non; T2m2_relax; T2m1_relax; T20_relax; T21_relax; T22_relax;   Tmn_non];
%         Tmn_timeSeries(3,:,:) = [T3m3_relax; T3m2_relax; T3m1_relax; T30_relax; T31_relax; T32_relax; T33_relax];
%         
%         Tmn_relax = Tmn_timeSeries(:,:,end);        
%         Tmn_relaxation = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift, obj.wShift_RMS, obj.wShift_FID, Tmn_relax);
%      end
     
    % % relaxation  with macroscopic anisotropy and chemical shift/B0 inhomogeneity (Cauchy/Lorentzian distribution for offset)
    % %             for Jen Model
    % % return: Tmn object after relaxation, 
    % %         Tmn Matrix with time evolution    
     function [Tmn_relaxation, Tmn_timeSeries] = relaxation_Jen(obj, tevo)
        Tmn_preRelax = obj.copy(); 
        Tmn_preRelax.Tmn = Tmn_preRelax.Tmn - Tmn_preRelax.Teq; % recover thermal equilibrium
 
        
        if tevo == 0
            ts = 0;
        else 
            ts = 0:tevo/obj.tsDataPoints:tevo;
        end
        Tmn_non = zeros([1 length(ts)]); % nonexisting Tmns
        T10_ts = ones([1 length(ts)]);  
  
        fm333  = getValues.get_f_Jen(-3, ts, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift_RMS,obj.w0);
        fm2kks = getValues.get_f_Jen(-2, ts, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift_RMS,obj.w0);
        fm1kks = getValues.get_f_Jen(-1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift_RMS,obj.w0);
        f0kks  = getValues.get_f_Jen( 0, ts, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift_RMS,obj.w0);
        f1kks  = getValues.get_f_Jen( 1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift_RMS,obj.w0);
        f2kks  = getValues.get_f_Jen( 2, ts, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift_RMS,obj.w0);
        f333   = getValues.get_f_Jen( 3, ts, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift_RMS,obj.w0);
          
%         disp(fm333(1,end))
        
        % zero quantum relaxation 
        f011 = f0kks(1,:); f013 = f0kks(2,:); f033 = f0kks(3,:); f022 = f0kks(4,:);
        T10_pre = Tmn_preRelax.getTmn(1, 0);
        T30_pre = Tmn_preRelax.getTmn(3, 0);
        T20_pre = Tmn_preRelax.getTmn(2, 0);
        % T10_30_relax = [f011 f013; f013; f033]*[T10_pre; T30_pre];
        % T10_relax = T10_30_relax(1);   T30_relax = T10_30_relax(2);
        T10_relax = f011*T10_pre+f013*T30_pre + T10_ts;   % add T10_ts to recover thermal equilibrium
        T30_relax = f013*T10_pre+f033*T30_pre;
        T20_relax = f022* T20_pre;
        
        % single quantum relaxation
        fm111 = fm1kks(1,:); fm113 = fm1kks(2,:); fm133 = fm1kks(3,:); fm122 = fm1kks(6,:);
        fm112 = -fm1kks(4,:); fm123 = -fm1kks(5,:);
        f111 = f1kks(1,:); f113 = f1kks(2,:); f133 = f1kks(3,:); f122 = f1kks(4,:);
        f112 = f1kks(4,:); f123 = f1kks(5,:);
        T1m1_pre = Tmn_preRelax.getTmn(1, -1); T11_pre = Tmn_preRelax.getTmn(1, 1);
        T3m1_pre = Tmn_preRelax.getTmn(3, -1); T31_pre = Tmn_preRelax.getTmn(3, 1);
        T2m1_pre = Tmn_preRelax.getTmn(2, -1); T21_pre = Tmn_preRelax.getTmn(2, 1);
        % n=-1
%         disp(size(exp(-1*(1i*obj.wShift)*ts)));
%         disp(size(cos(2*pi*obj.wShift_Freq.*ts)));
        T1m1_relax = (fm111*T1m1_pre + fm112*T2m1_pre + fm113*T3m1_pre).*exp(-1*(1i*obj.wShift)*ts);   %.*cos(2*pi*obj.wShift_Freq*ts)
        T3m1_relax = (fm113*T1m1_pre + fm123*T2m1_pre + fm133*T3m1_pre).*exp(-1*(1i*obj.wShift)*ts);
        T2m1_relax = (fm112*T1m1_pre + fm122*T2m1_pre + fm123*T3m1_pre).*exp(-1*(1i*obj.wShift)*ts);
        % n=+1
        T11_relax = (f111*T11_pre + f112*T21_pre + f113*T31_pre).*exp(1*(1i*obj.wShift)*ts);   
        T31_relax = (f113*T11_pre + f123*T21_pre + f133*T31_pre).*exp(1*(1i*obj.wShift)*ts);
        T21_relax = (f112*T11_pre + f122*T21_pre + f123*T31_pre).*exp(1*(1i*obj.wShift)*ts);
        
        % multiple quantum relaxation
        % m = 2
        f222 = f2kks(1,:); f233 = f2kks(2,:); f223 = f2kks(3,:); 
        fm222 = fm2kks(1,:); fm233 = fm2kks(2,:); fm223 = f2kks(3,:);
        T2m2_pre = Tmn_preRelax.getTmn(2, -2); T22_pre = Tmn_preRelax.getTmn(2, 2);
        T3m2_pre = Tmn_preRelax.getTmn(3, -2); T32_pre = Tmn_preRelax.getTmn(3, 2);
        % n=-2
        T2m2_relax = (fm222*T2m2_pre + fm223*T3m2_pre).*exp(-2*(1i*obj.wShift)*ts);
        T3m2_relax = (fm233*T3m2_pre + fm223*T2m2_pre).*exp(-2*(1i*obj.wShift)*ts);
        % n=+2
        T22_relax = (f222*T22_pre + f223*T32_pre).*exp(2*(1i*obj.wShift)*ts);
        T32_relax = (f233*T32_pre + f223*T22_pre).*exp(2*(1i*obj.wShift)*ts);
        % m = 3 
        T3m3_pre = Tmn_preRelax.getTmn(3, -3); T33_pre = Tmn_preRelax.getTmn(3, 3);
        T3m3_relax = fm333* T3m3_pre.*exp(-3*(1i*obj.wShift)*ts);
        T33_relax  = f333* T33_pre.*exp(3*(1i*obj.wShift)*ts);

        Tmn_timeSeries = zeros([3 7 length(ts)]);
        Tmn_timeSeries(1,:,:) = [   Tmn_non;    Tmn_non; T1m1_relax; T10_relax; T11_relax;   Tmn_non;   Tmn_non];
        Tmn_timeSeries(2,:,:) = [   Tmn_non; T2m2_relax; T2m1_relax; T20_relax; T21_relax; T22_relax;   Tmn_non];
        Tmn_timeSeries(3,:,:) = [T3m3_relax; T3m2_relax; T3m1_relax; T30_relax; T31_relax; T32_relax; T33_relax];
        
        Tmn_relax = Tmn_timeSeries(:,:,end);        
        Tmn_relaxation = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift, obj.wShift_RMS, obj.wShift_FID, Tmn_relax);
     end
     
     
    % % return FID --> only Relaxation in Tmn(1,-1) channel
    % % receiver phase phaseRX     
    function [FID] = getFIDphase(obj, ts, phaseRX)
        Tmn_preRelax = obj.Tmn; % für T10 Relax muss Teq abgezogen und am Ende wieder addiert werden
        fm1kks = getValues.get_f_Jen(-1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift_RMS, obj.w0);
        % f1kks = getValues.get_f(1, ts, obj.tauC, obj.wQ, obj.w0);
        % single quantum relaxation
        fm111 = fm1kks(1,:); fm113 = fm1kks(2,:); fm112 = -fm1kks(4,:);
        % f111 = f1kks(1); f113 = f1kks(2); 
        T1m1_pre = TmnEvo.getTmn_matrix(Tmn_preRelax, 1, -1); % T11_pre = Tmn_preRelax(1, 1);
        T3m1_pre = TmnEvo.getTmn_matrix(Tmn_preRelax, 3, -1); % T31_pre = Tmn_preRelax(1, 1);
        T2m1_pre = TmnEvo.getTmn_matrix(Tmn_preRelax, 2, -1);
        
%         T1m1_relax = fm111*T1m1_pre+fm113*T3m1_pre;   
        % T11_relax = f111*T11_pre+f113*T31_pre;
        
%         T1m1_relax = (fm111*T1m1_pre + fm112*T2m1_pre + fm113*T3m1_pre) .* exp(-1*(obj.wShift_FID + 1i*obj.wShift)*ts); 
        T1m1_relax = (fm111*T1m1_pre + fm112*T2m1_pre + fm113*T3m1_pre); 
        
        FID = T1m1_relax * exp(1i *deg2rad(phaseRX));

%         figure(); plot(ts, real(FID)); 
    end


    function [FID,times] = relaxationFID(obj,times)
        tsDataPoints   = 2048
        [Tmn_evo, Tmn_timeSeries_evo_1] = obj.relaxation_Jen(times);
        T11relax = Tmn_timeSeries_evo_1(1,5,:);
        
%         for idx = 1:length(times)
%             ts = times(idx);
%             fm1kks = getValues.get_f_Jen(-1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift_RMS, obj.w0);            
%         end
    end
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     
     function [Tmn_softPulse, Tmn_timeSeries, ts] = softPulse(obj, t_duration, pulseFun, pulsePhase)
         tsim           = linspace(0.0, t_duration, obj.noFrames);
         
         MStart         = TmnEvo.Tmn2M(obj.Tmn);
         M0             = TmnEvo.Tmn2M( TmnEvo.Teq);
         
         J0             = getValues.get_J_Jen(0, obj.tauC, obj.wQ, obj.Jen, obj.w0);
         J1             = getValues.get_J_Jen(1, obj.tauC, obj.wQ, obj.Jen, obj.w0);
         J2             = getValues.get_J_Jen(2, obj.tauC, obj.wQ, obj.Jen, obj.w0);
                  
         qVec           = [0 1 1 0 1 1 2 2 0 1 1 2 2 3 3 ];
         M_B0_inhomo    = diag(-qVec.*obj.wShift_RMS);
         
         sodium         = @(J0, J1, J2, Delta, omegaQ, omega1, phi, M, M0, t) ([(-2/5).*J1+(-8/5).*J2,omega1(t).*sin(-phi(t)),(1i*(-1)).*omega1(t).*cos(phi(t)),0,0,0,0,0,(-4/5).*J1+(4/5).*J2,0,0,0,0,0,0;
                (-1).*omega1(t).*sin(-phi(t)),(-3/5).*J0+(-1).*J1+(-2/5).*J2,Delta*1i,0,0,1i.*(3/5).^(1/2).*omegaQ,0,0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,0,0;
                (1i*(-1)).*omega1(t).*cos(phi(t)),Delta*1i,(-3/5).*J0+(-1).*J1+(-2/5).*J2,0,1i.*(3/5).^(1/2).*omegaQ,0,0,0,0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,0;
                0,0,0,(-2).*J1+(-2).*J2,3.^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*3.^(1/2).*omega1(t).*cos(phi(t)),0,0,0,0,0,0,0,0,0;
                0,0,1i.*(3/5).^(1/2).*omegaQ,(-1).*3.^(1/2).*omega1(t).*sin(-phi(t)),(-1).*J0+(-1).*J1+(-2).*J2,Delta*1i,(1i*(-1)).*omega1(t).*cos(phi(t)),omega1(t).*sin(-phi(t)),0,0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,0;
                0,1i.*(3/5).^(1/2).*omegaQ,0,(1i*(-1)).*3.^(1/2).*omega1(t).*cos(phi(t)),Delta*1i,(-1).*J0+(-1).*J1+(-2).*J2,omega1(t).*sin(-phi(t)),(1i*(-1)).*omega1(t).*cos(phi(t)),0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,0,0;
                0,0,0,0,(1i*(-1)).*omega1(t).*cos(phi(t)),(-1).*omega1(t).*sin(-phi(t)),(-1).*J0+(-2).*J1+(-1).*J2,2*Delta*1i,0,0,0,0,1i.*omegaQ,0,0;
                0,0,0,0,(-1).*omega1(t).*sin(-phi(t)),(1i*(-1)).*omega1(t).*cos(phi(t)),2*Delta*1i,(-1).*J0+(-2).*J1+(-1).*J2,0,0,0,1i.*omegaQ,0,0,0;
                (-4/5).*J1+(4/5).*J2,0,0,0,0,0,0,0,(-8/5).*J1+(-2/5).*J2,6.^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*6.^(1/2).*omega1(t).*cos(phi(t)),0,0,0,0;
                0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,1i.*(2/5).^(1/2).*omegaQ,0,0,(-1).*6.^(1/2).*omega1(t).*sin(-phi(t)),(-2/5).*J0+(-1).*J1+(-3/5).*J2,Delta*1i,(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),(5/2).^(1/2).*omega1(t).*sin(-phi(t)),0,0;
                0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,(1i*(-1)).*6.^(1/2).*omega1(t).*cos(phi(t)),Delta*1i,(-2/5).*J0+(-1).*J1+(-3/5).*J2,(5/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),0,0;
                0,0,0,0,0,0,0,1i.*omegaQ,0,(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),(-1).*(5/2).^(1/2).*omega1(t).*sin(-phi(t)),(-1).*J0+(-1).*J2,2*Delta*1i,(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t)),(3/2).^(1/2).*omega1(t).*sin(-phi(t));
                0,0,0,0,0,0,1i.*omegaQ,0,0,(-1).*(5/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),2*Delta*1i,(-1).*J0+(-1).*J2,(3/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t));
                0,0,0,0,0,0,0,0,0,0,0,(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t)),(-1).*(3/2).^(1/2).*omega1(t).*sin(-phi(t)),(-1).*J1+(-1).*J2,3*Delta*1i;
                0,0,0,0,0,0,0,0,0,0,0,(-1).*(3/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t)),3*Delta*1i,(-1).*J1+(-1).*J2] + M_B0_inhomo) * M...
                + [2/5.*J1+8/5.*J2;0;0;0;0;0;0;0;4/5.*J1-4/5.*J2;0;0;0;0;0;0];

        % numerical integration of the Equations
        options             = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,15));    % solver accuracy
        [ts, M_softPulse_timeSeries] = ode45(@(t, M)sodium(J0, J1, J2, 2*pi*obj.DeltaHz, 2*pi*obj.wQbar, pulseFun, pulsePhase, M, M0, t_duration), tsim, MStart, options); 

        Tmn_timeSeries      = TmnEvo.Ms2Tmns(M_softPulse_timeSeries);
        M_softPulse         = M_softPulse_timeSeries(end, :);
        
        Tmn_softPulseMatrix = TmnEvo.M2Tmn(M_softPulse);
        Tmn_softPulse       = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, obj.wShift, obj.wShift_RMS, obj.wShift_FID, Tmn_softPulseMatrix); 
     end
     
     function [Tmn_softPulse, Tmn_timeSeries, ts] = softPulse_rect(obj, t_duration, pulseStrength, pPhase)
         pulse           = @(t)  0.0 + (t >= 0) .* (t <= t_duration) .* pulseStrength;            
         tr_phase        = @(t)  0.0 + (t >= 0) .* (t <= t_duration) .* deg2rad(pPhase);         
         [Tmn_softPulse, Tmn_timeSeries, ts] = softPulse(obj, t_duration, pulse, tr_phase);
     end
     
     
     
     
     
     
     
     
end
end
