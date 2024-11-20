classdef PhaseCyclesDZ < handle
    properties
        B0 = 9.4 % T
        
        w0
        tauC
        wQ
        wQbar          = 0;
        Jen
        
        wShift         = 0;         % B-field inhomogeneity -> mean(Omega)
        wShift_RMS     = 0;         % B-field inhomogeneity -> ~std(Omega)  --> intrinsical, no refocus by 180
        wShift_dist
        FreqShift      = 0;
        
        flipShift_dist90
        wShift_FID     = 0;         % B-field inhomogeneity for 1. Dim FID  --> refocus by 180 possible (but not performed in FID
            
        tevoStep       = 0.1000 * 1e-3; % s
        tevo0          = 0.1305 * 1e-3; % s 
        tmix           = 0.1500 * 1e-3; % s
        tmixStep       = 0.2 * 1e-3; % s
    
        startPhase     = 90; % degrees
        alphaStep      = 45 %degrees
        alphas 
    
        nSpins         = 100;
    
        dur1           = 50e-6; %dur1 as defined in TQTTPI .ppg file  
    
        NumPhaseCycles 
        PulseAngles    = []
        
        TR             = 300e-3; % s  repetition time
        
    % %     timeFID        = 200e-3; % s
        dataPoints     = 2048;
        
        deadtimeFID    =   6e-6; % s
        dwelltimeFID   =  25e-6; % s
    
        noisePower     = 0; % noise Power in Watts, linear, currently only for T3i sequence
        
        
        flip90         = 90; % can be set to different values to simulate B1+ 
        flip180        = 180;% inhomogeneities
        
        TEvec          = 0;
        TE0            = 0;
        NumTE          = 256;  
        TEstep         = 0.1e-3;
        TEmax          = 400e-3;
        
        TmnVisMethod   = 0;
        TmnVisPulse    = false;
        
        second180phase = 0; %only important for
    end
methods(Static)  
    % % general 3 Pulse sequence phase step
    % % phase of 3rd pulse is always 0 
    function [TmnEnd, TQ, SQ, Tmn_timeSeries,Tmn_timeSeries_relax] = general3Pulse( TmnStart, alpha, beta, tevo, tmix, TR, flip1, flip2, flip3,phaseRX)
        [Tmn_pulse90_1] = TmnStart.pulse(flip1, alpha);
        [Tmn_evo, Tmn_timeSeries_evo_1] = Tmn_pulse90_1.relaxation_Jen(tevo);
%         Tmn_evo.setTmn(3,  1, 0);
%         Tmn_evo.setTmn(3, -1, 0);
%         Tmn_evo.setTmn(1,  1, 0);
%         Tmn_evo.setTmn(1, -1, 0);
%         Tmn_evo.setTmn(3,  0, 0);
%         Tmn_evo.setTmn(1,  0, 0);   
        [Tmn_pulse90_2] = Tmn_evo.pulse(flip2, beta);    
%         Tmn_pulse90_2.setTmn(3,  1, 0);
%         Tmn_pulse90_2.setTmn(3, -1, 0);
%         Tmn_pulse90_2.setTmn(3,  3, 0);
%         Tmn_pulse90_2.setTmn(3, -3, 0);
%         Tmn_pulse90_2.setTmn(1,  1, 0);
%         Tmn_pulse90_2.setTmn(1, -1, 0);
%         Tmn_pulse90_2.setTmn(3,  0, 0);
%         Tmn_pulse90_2.setTmn(1,  0, 0);   
        [Tmn_evo_2, Tmn_timeSeries_mix] = Tmn_pulse90_2.relaxation_Jen(tmix);
%         Tmn_evo_2.setTmn(3,  3, 0);
%         Tmn_evo_2.setTmn(3, -3, 0);
%         Tmn_evo_2.setTmn(3,  1, 0);
%         Tmn_evo_2.setTmn(3, -1, 0);
        [Tmn_pulse90_3] = Tmn_evo_2.pulse(flip3, 0);
%         Tmn_pulse90_3.setTmn(1, -1, 0);
%         Tmn_pulse90_3.setTmn(1, 1, 0)
%         Tmn_pulse90_2.setTmn(1,  0, 0);

        TmnEnd = Tmn_pulse90_3.copy();
%         TmnEnd = Tmn_pulse90_2;
%         TmnEnd = Tmn_evo;

        %create complete timeseries matrix
        [Tmn_evo_4, Tmn_timeSeries_relax] = Tmn_pulse90_3.relaxation_Jen(TR); % TR decides how far the relaxation is sampled
        Tmn_timeSeries_relax = Tmn_timeSeries_relax * exp(1i *deg2rad(phaseRX));

        Tmn_timeSeries = cat(3,Tmn_timeSeries_evo_1(:,:,1:end-1), Tmn_timeSeries_mix(:,:,1:end-1), Tmn_timeSeries_relax(:,:,1:end-1));
        
        % calculate TQ signal (T3-3, T33 after mixing period)
        TQp = Tmn_evo_2.getTmn(3, +3)*getValues.get_WignerD(3, -1, +3, flip3); 
        TQm = Tmn_evo_2.getTmn(3, -3)*getValues.get_WignerD(3, -1, -3, flip3);
        TQ = TQp+TQm;
        
        % calculate SQ signal (T1-1, T11, T3-1, T31 after mixing period)
        SQp = Tmn_evo_2.getTmn(1, +1)*getValues.get_WignerD(1, -1, +1, flip3) + Tmn_evo_2.getTmn(3, +1)*getValues.get_WignerD(3, -1, +1, flip3); 
        SQm = Tmn_evo_2.getTmn(1, -1)*getValues.get_WignerD(1, -1, -1, flip3) + Tmn_evo_2.getTmn(3, -1)*getValues.get_WignerD(3, -1, -1, flip3); 
        SQ = SQp + SQm;
    end  
    
    function [TmnEnd, TQ, SQ, Tmn_timeSeries,Tmn_timeSeries_relax, Tmn_timeSeries_evo1ADC] = general3PulseWithADC1( TmnStart, alpha, beta, tevo, tmix, TR, flip1, flip2, flip3,phaseRX)
        [Tmn_pulse90_1] = TmnStart.pulse(flip1, alpha);
        tevoAdapted = tevo - 10e-06;
        % tevoAdapted = tevo - obj.dur1 - 10e-06;
        [Tmn_evo, Tmn_timeSeries_evo_1] = Tmn_pulse90_1.relaxation_Jen(tevo);
        [Tmn_evoADC1, Tmn_timeSeries_evo1ADC] = Tmn_pulse90_1.relaxation_Jen(tevoAdapted);

%         Tmn_evo.setTmn(3,  1, 0);
%         Tmn_evo.setTmn(3, -1, 0);
%         Tmn_evo.setTmn(1,  1, 0);
%         Tmn_evo.setTmn(1, -1, 0);
%         Tmn_evo.setTmn(3,  0, 0);
%         Tmn_evo.setTmn(1,  0, 0);   
        [Tmn_pulse90_2] = Tmn_evo.pulse(flip2, beta);    
%         Tmn_pulse90_2.setTmn(3,  1, 0);
%         Tmn_pulse90_2.setTmn(3, -1, 0);
%         Tmn_pulse90_2.setTmn(3,  3, 0);
%         Tmn_pulse90_2.setTmn(3, -3, 0);
%         Tmn_pulse90_2.setTmn(1,  1, 0);
%         Tmn_pulse90_2.setTmn(1, -1, 0);
%         Tmn_pulse90_2.setTmn(3,  0, 0);
%         Tmn_pulse90_2.setTmn(1,  0, 0);   
        [Tmn_evo_2, Tmn_timeSeries_mix] = Tmn_pulse90_2.relaxation_Jen(tmix);
%         Tmn_evo_2.setTmn(3,  3, 0);
%         Tmn_evo_2.setTmn(3, -3, 0);
%         Tmn_evo_2.setTmn(3,  1, 0);
%         Tmn_evo_2.setTmn(3, -1, 0);
        [Tmn_pulse90_3] = Tmn_evo_2.pulse(flip3, 0);
%         Tmn_pulse90_3.setTmn(1, -1, 0);
%         Tmn_pulse90_3.setTmn(1, 1, 0)
%         Tmn_pulse90_2.setTmn(1,  0, 0);

        TmnEnd = Tmn_pulse90_3.copy();
%         TmnEnd = Tmn_pulse90_2;
%         TmnEnd = Tmn_evo;

        %create complete timeseries matrix
        [Tmn_evo_4, Tmn_timeSeries_relax] = Tmn_pulse90_3.relaxation_Jen(TR); % TR decides how far the relaxation is sampled
        Tmn_timeSeries_relax = Tmn_timeSeries_relax * exp(1i *deg2rad(phaseRX));

        Tmn_timeSeries = cat(3,Tmn_timeSeries_evo_1(:,:,1:end-1), Tmn_timeSeries_mix(:,:,1:end-1), Tmn_timeSeries_relax(:,:,1:end-1));
        
        % calculate TQ signal (T3-3, T33 after mixing period)
        TQp = Tmn_evo_2.getTmn(3, +3)*getValues.get_WignerD(3, -1, +3, flip3); 
        TQm = Tmn_evo_2.getTmn(3, -3)*getValues.get_WignerD(3, -1, -3, flip3);
        TQ = TQp+TQm;
        
        % calculate SQ signal (T1-1, T11, T3-1, T31 after mixing period)
        SQp = Tmn_evo_2.getTmn(1, +1)*getValues.get_WignerD(1, -1, +1, flip3) + Tmn_evo_2.getTmn(3, +1)*getValues.get_WignerD(3, -1, +1, flip3); 
        SQm = Tmn_evo_2.getTmn(1, -1)*getValues.get_WignerD(1, -1, -1, flip3) + Tmn_evo_2.getTmn(3, -1)*getValues.get_WignerD(3, -1, -1, flip3); 
        SQ = SQp + SQm;
    end  



    % % general 3 Pulse sequence phase step
    % % phase of 3rd pulse is always 0 
    function [TmnEnd, TQ, SQ] = general3Pulse3phase( TmnStart, alpha, beta, phase3, tevo, tmix, flip1, flip2, flip3)
        [Tmn_pulse90_1] = TmnStart.pulse(flip1, alpha);
        [Tmn_evo, Tmn_timeSeries_evo_1] = Tmn_pulse90_1.relaxation_Jen(tevo);
        [Tmn_pulse90_2] = Tmn_evo.pulse(flip2, beta);      
        [Tmn_evo_2, Tmn_timeSeries_mix] = Tmn_pulse90_2.relaxation_Jen(tmix);
        [Tmn_pulse90_3] = Tmn_evo_2.pulse(flip3, phase3);

        TmnEnd = Tmn_pulse90_3;

        % calculate TQ signal (T3-3, T33 after mixing period)
        TQp = Tmn_evo_2.getTmn(3, +3)*getValues.get_WignerD(3, -1, +3, flip3); 
        TQm = Tmn_evo_2.getTmn(3, -3)*getValues.get_WignerD(3, -1, -3, flip3);
        TQ = TQp+TQm;
        
        % calculate SQ signal (T1-1, T11, T3-1, T31 after mixing period)
        SQp = Tmn_evo_2.getTmn(1, +1)*getValues.get_WignerD(1, -1, +1, flip3) + Tmn_evo_2.getTmn(3, +1)*getValues.get_WignerD(3, -1, +1, flip3); 
        SQm = Tmn_evo_2.getTmn(1, -1)*getValues.get_WignerD(1, -1, -1, flip3) + Tmn_evo_2.getTmn(3, -1)*getValues.get_WignerD(3, -1, -1, flip3); 
        SQ = SQp + SQm;
    end  
    
    
    
    % % general 3 Pulse sequence phase step with first pulse being a
    % % rectagular soft pulse of length tevo
    % % phase of 3rd pulse is always 0 
    % % 
    function [TmnEnd, TQ, SQ] = general3Pulse_softPulse( TmnStart, alpha, beta, tevo, tmix, pulseStrength, flip2, flip3)
        [Tmn_softPulse, Tmn_timeSeries_softPulse, ts] = TmnStart.softPulse_rect(tevo, pulseStrength, alpha);
        [Tmn_pulse90_2] = Tmn_softPulse.pulse(flip2, beta);
        [Tmn_evo_2, Tmn_timeSeries_mix] = Tmn_pulse90_2.relaxation_Jen(tmix);
        [Tmn_pulse90_3] = Tmn_evo_2.pulse(flip3, 0);

        TmnEnd = Tmn_pulse90_3;
        
        % calculate TQ signal (T3-3, T33 after mixing period)
        TQp = Tmn_evo_2.getTmn(3, +3)*getValues.get_WignerD(3, -1, +3, flip3); 
        TQm = Tmn_evo_2.getTmn(3, -3)*getValues.get_WignerD(3, -1, -3, flip3);
        TQ = TQp+TQm;
        
        % calculate SQ signal (T1-1, T11, T3-1, T31 after mixing period)
        SQp = Tmn_evo_2.getTmn(1, +1)*getValues.get_WignerD(1, -1, +1, flip3) + Tmn_evo_2.getTmn(3, +1)*getValues.get_WignerD(3, -1, +1, flip3); 
        SQm = Tmn_evo_2.getTmn(1, -1)*getValues.get_WignerD(1, -1, -1, flip3) + Tmn_evo_2.getTmn(3, -1)*getValues.get_WignerD(3, -1, -1, flip3); 
        SQ = SQp + SQm;
    end 
    
    % % general 3 Pulse sequence phase step with 180° refocussing pulse
    % % phase of 3rd pulse is always 0 
    function [TmnEnd, TQ, SQ, Tmn_timeSeries,Tmn_timeSeries_relax] = general3Pulse_w180(TmnStart, alpha, beta, tevo, tmix, TR, flip1, flip2, flip3, flipRefocus, phaseRX)
        [Tmn_pulse90_1] = TmnStart.pulse(flip1, alpha);
        [Tmn_evo_1, Tmn_timeSeries_evo_1] = Tmn_pulse90_1.relaxation_Jen(tevo/2);
        [Tmn_pulse180] = Tmn_evo_1.pulse(flipRefocus, alpha+90); %normal: alpha+90
        [Tmn_evo_2, Tmn_timeSeries_evo_2] = Tmn_pulse180.relaxation_Jen(tevo/2);
        [Tmn_pulse90_2] = Tmn_evo_2.pulse(flip2, beta);
        [Tmn_evo_3, Tmn_timeSeries_mix] = Tmn_pulse90_2.relaxation_Jen(tmix);
        [Tmn_pulse90_3] = Tmn_evo_3.pulse(flip3, 0);
        TmnEnd = Tmn_pulse90_3.copy();
        
        %create complete timeseries matrix
        [~, Tmn_timeSeries_relax] =Tmn_pulse90_3.relaxation_Jen(TR);
        Tmn_timeSeries_relax = Tmn_timeSeries_relax * exp(1i *deg2rad(phaseRX));
%         disp(size(Tmn_timeSeries_relax))
        Tmn_timeSeries = cat(3,Tmn_timeSeries_evo_1, Tmn_timeSeries_evo_2, Tmn_timeSeries_mix, Tmn_timeSeries_relax);
%         disp(size(Tmn_timeSeries_evo_1))

        % calculate TQ signal (T3-3, T33 after mixing period)
        TQp = Tmn_evo_3.getTmn(3, +3)*getValues.get_WignerD(3, -1, +3, flip3); 
        TQm = Tmn_evo_3.getTmn(3, -3)*getValues.get_WignerD(3, -1, -3, flip3);
        TQ = TQp+TQm;
        
        % calculate SQ signal (T1-1, T11, T3-1, T31 after mixing period)
        SQp = Tmn_evo_3.getTmn(1, +1)*getValues.get_WignerD(1, -1, +1, flip3) + Tmn_evo_3.getTmn(3, +1)*getValues.get_WignerD(3, -1, +1, flip3); 
        SQm = Tmn_evo_3.getTmn(1, -1)*getValues.get_WignerD(1, -1, -1, flip3) + Tmn_evo_3.getTmn(3, -1)*getValues.get_WignerD(3, -1, -1, flip3); 
        SQ = SQp + SQm;
    end
    
    % % general 3 Pulse sequence with 180° between 2nd and 3rd pulse
    % Dominik
    function [TmnEnd, TQ, SQ, Tmn_timeSeries,Tmn_timeSeries_relax] = general3Pulse_w180in23(TmnStart, alpha, beta, tevo, tmix, TR, flip1, flip2, flip3, flipRefocus, phaseRX)
        [Tmn_pulse90_1] = TmnStart.pulse(flip1, alpha);
        [Tmn_evo_1, Tmn_timeSeries_evo_1] = Tmn_pulse90_1.relaxation_Jen(tevo);
        
        [Tmn_pulse90_2] = Tmn_evo_1.pulse(flip2, beta);

        [Tmn_mix_1, Tmn_timeSeries_mix1] = Tmn_pulse90_2.relaxation_Jen(tmix/2);
        [Tmn_pulse180] = Tmn_mix_1.pulse(flipRefocus, beta+180); %normal: alpha+90
        [Tmn_mix_2, Tmn_timeSeries_mix2] = Tmn_pulse180.relaxation_Jen(tmix/2);
        [Tmn_pulse90_3] = Tmn_mix_2.pulse(flip3, 0);
        TmnEnd = Tmn_pulse90_3.copy();
        
        %create complete timeseries matrix
        [~, Tmn_timeSeries_relax] =Tmn_pulse90_3.relaxation_Jen(TR);
        Tmn_timeSeries_relax = Tmn_timeSeries_relax * exp(1i *deg2rad(phaseRX));
%         disp(size(Tmn_timeSeries_relax))
        Tmn_timeSeries = cat(3,Tmn_timeSeries_evo_1, Tmn_timeSeries_mix1, Tmn_timeSeries_mix2, Tmn_timeSeries_relax);
%         disp(size(Tmn_timeSeries_evo_1))

        % calculate TQ signal (T3-3, T33 after mixing period)
        TQp = Tmn_mix_2.getTmn(3, +3)*getValues.get_WignerD(3, -1, +3, flip3); 
        TQm = Tmn_mix_2.getTmn(3, -3)*getValues.get_WignerD(3, -1, -3, flip3);
        TQ = TQp+TQm;
        
        % calculate SQ signal (T1-1, T11, T3-1, T31 after mixing period)
        SQp = Tmn_mix_2.getTmn(1, +1)*getValues.get_WignerD(1, -1, +1, flip3) + Tmn_mix_2.getTmn(3, +1)*getValues.get_WignerD(3, -1, +1, flip3); 
        SQm = Tmn_mix_2.getTmn(1, -1)*getValues.get_WignerD(1, -1, -1, flip3) + Tmn_mix_2.getTmn(3, -1)*getValues.get_WignerD(3, -1, -1, flip3); 
        SQ = SQp + SQm;
    end
    
    % % correct phase of 1st dimenstion FIDs
    function [FID] = phase_cor_FID(fid)
%         spectrum = fftshift(fft(fid(1,:)));
%         [~,pos] = max(abs(spectrum));
%         banana = @(x) -real(spectrum(pos)*exp(1i*x));
%         [x(1),~] = fminsearch(banana,[0]);
%         if max(real(fid(1,:)* exp(1i*x(1)))) < 0
%             FID = fid * exp(-1i*x(1));
%         else
%             FID = fid * exp(1i*x(1));
%         end

        % % alternative phase correction method without fft 
        fidi = fid(1,:);
        banana = @(x) -mean(real(fidi*exp(1i*x)));
        [x(1),~] = fminsearch(banana,[0]);
        if max(real(fid(1,:)* exp(1i*x(1)))) < 0
            FID = fid * exp(-1i*x(1));
        else
            FID = fid * exp(1i*x(1));
        end

% %         % % alternative phase correction: 
%         fidi = mean(fid,2);
%         banana = @(x) mean(abs(imag(fidi*exp(1i*x))));
%         [x(1),~] = fminsearch(banana,[0]);
%         if max(real(fid(1,:)* exp(1i*x(1)))) < 0
%             FID = fid * exp(-1i*x(1));
%         else 
%             FID = fid * exp(1i*x(1));
%         end
    end
    
    % % TQTPPI FID
    % % both FIDs are added for DQ suppression
    % % Peak hight of 1st dimension spectra used for 2nd dimension FID
    function [tevos, TQTPPI_FID] = get_TQTPPI_FID(tevos, FIDs_p1, FIDs_p2)
        sizeP1 = size(FIDs_p1);
        FIDs = FIDs_p1 + FIDs_p2;

        % DC correction: mean of all 1.dim FIDs last data point should be = 0
%         FIDs = FIDs - mean(FIDs(:,end));   
        
%         FIDs = PhaseCycle.phase_cor_FID(FIDs);
        
%         FIDs_ft =  zeros(sizeP1);
%         for j = 1:sizeP1(1)
%             FID = FIDs(j,:);
%             spec = fftshift(fft(FID.'));
%             FIDs_ft(j,:) =   spec;
% %             FIDs_ft(j,:) =  abs( spec);
% %             FIDs_ft(j,:) =   FID.';
%         end
%         maxFIDs = FIDs_ft(:,fix(sizeP1(2)/2)+1);
%         TQTPPI_FID = maxFIDs; 
        disp(size(FIDs));
        sumFIDs = sum(FIDs,2);
        TQTPPI_FID = sumFIDs;
 
%         TQTPPI_FID = TQTPPI_FID - mean(TQTPPI_FID);  
    end 
 
    % % TQTPPI FID
    % % both FIDs are added for DQ suppression
    % % Peak hight of 1st dimension spectra used for 2nd dimension FID
    % % filter out contributions from exp-exp biexponential Terms
    function [tevos, TQTPPI_FID] = get_TQTPPI_FID_woTQ(tevos, FIDs_p1, FIDs_p2)
        sizeP1 = size(FIDs_p1);
        FIDs = FIDs_p1 + FIDs_p2;

        % DC correction: mean of all 1.dim FIDs last data point should be = 0
%         FIDs = FIDs - mean(FIDs(:,end));   
        
%         FIDs = PhaseCycle.phase_cor_FID(FIDs);
        
        FIDs_ft =  zeros(sizeP1);
        for j = 1:sizeP1(1)
            FID = FIDs(j,:);
            spec = fftshift(fft(FID.'));
            FIDs_ft(j,:) =   spec;
%             FIDs_ft(j,:) =   spec;
%             FIDs_ft(j,:) =   FID.';
        end
%         maxFIDs = FIDs_ft(:,fix(sizeP1(2)/2)+1);
%         TQTPPI_FID = maxFIDs; 
        
        sumFIDs = sum(FIDs_ft,2);
        TQTPPI_FID = sumFIDs;
 
%         TQTPPI_FID = TQTPPI_FID - mean(TQTPPI_FID);  
    end    


    
end
methods
    % constructor
    function obj = PhaseCyclesDZ(B0, tauC, wQ, wQbar, Jen, wShift, wShift_RMS, wShift_FID )
        obj.wQbar          = wQbar;
        obj.B0             = B0;
        obj.tauC           = tauC;
        obj.wQ             = wQ;
        obj.alphas         = (0:obj.alphaStep:359) + obj.startPhase; % degrees   
        obj.Jen            = Jen;
        
        obj.wShift         = wShift;
        obj.wShift_RMS     = wShift_RMS;
        obj.NumPhaseCycles = 100;
        
        obj.wShift_FID     = wShift_FID;

    end   

    
    % % 90 degree pulse with wShift only relax
    function [ FID, times] = Pulse90hard_f111(obj)
        times = obj.deadtimeFID + obj.dwelltimeFID.*(0:obj.dataPoints-1);

        wShift = 0;
        wShift_RMS = 0;
        wShiftFID = 0;
        Tmn_i = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, wShift, wShift_RMS, wShiftFID);
        [Tmn_pulse1] = Tmn_i.pulse(obj.flip90, 0);
        f1kks  = getValues.get_f_Jen( 1, times, Tmn_pulse1.tauC, Tmn_pulse1.wQ, Tmn_pulse1.wQbar, Tmn_pulse1.Jen, Tmn_pulse1.wShift_RMS,Tmn_pulse1.w0);
%         size(f1kks)
        
        wShift_relax = random(obj.wShift_dist);
%         wShift_relax = randn(1)*4;
        FID = f1kks(1,:).*exp(1i*wShift_relax*times);
        
    end

    %%%%%%%% Sequences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fixed TQTPPI wo. 180 
    function [ times,FID,FIDs] = TQTPPIfixed_wo180_VJ(obj)
        fprintf('Here! \n');
        times = obj.deadtimeFID +(0:obj.dataPoints-1) .* (obj.TR/obj.dataPoints);
        %times = times .* 2;
        lengthAlpha = length(obj.alphas);
        tevo = obj.tevo0;
        nFIDs = obj.NumPhaseCycles * lengthAlpha;
        
        FIDs = zeros(nFIDs,obj.nSpins,length(times),2);
        for j = 1:obj.NumPhaseCycles
            for alphaIndex = 1:lengthAlpha
                idx = (j-1)*lengthAlpha + alphaIndex ;
                alpha = obj.alphas(alphaIndex);
                beta_p1 = alpha + 90;
                beta_p2 = alpha - 90;
%                 beta_p1 = alpha + 0;
%                 beta_p2 = alpha + 180;
        
                
                FIDs_p1 = zeros(obj.nSpins,length(times));
                FIDs_p2 = zeros(obj.nSpins,length(times));
                %fprintf('AlphaIdx: %d \n', alphaIndex);
                if obj.nSpins > 50

                    parfor idxSpin = 1:obj.nSpins 
                        w_Shift = obj.FreqShift + random(obj.wShift_dist);
  %                     w_Shift = 0;
    
                        TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                        [~,~,~,~,Tmn_p1_Spin_relax] = PhaseCyclesDZ.general3Pulse(TmnStart, alpha, beta_p1, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,0);
                        [~,~,~,~,Tmn_p2_Spin_relax] = PhaseCyclesDZ.general3Pulse(TmnStart, alpha, beta_p2, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,0);
    
                        FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                        FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                    end
                else
                    for idxSpin = 1:obj.nSpins 
                        w_Shift = obj.FreqShift + random(obj.wShift_dist);
  %                     w_Shift = 0;
    
                        TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                        [~,~,~,~,Tmn_p1_Spin_relax] = PhaseCyclesDZ.general3Pulse(TmnStart, alpha, beta_p1, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,0);
                        [~,~,~,~,Tmn_p2_Spin_relax] = PhaseCyclesDZ.general3Pulse(TmnStart, alpha, beta_p2, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,0);
    
                        FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                        FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                    end
                end
                    
                
                FIDs(idx,:,:,1) = FIDs_p1;
                FIDs(idx,:,:,2) = FIDs_p2;
                
                if mod(idx, 10) == 0
                    fprintf('PhaceCycle %d/%d completed\n',idx,nFIDs);
                end
            end
        end
                   
        FID_sumSpins = sum(FIDs,2);
        FID_sumFIDs = sum(FID_sumSpins,3);
        FID_sumBetas = FID_sumFIDs(:,:,:,1)+FID_sumFIDs(:,:,:,2);
%         FID_sumBetas = sum(FID_sumFIDs,4);
%         FID = reshape(FID_sumBetas,nFIDs,1);
        
        FID = getValues.switchReIm(FID_sumBetas);
    end
    % TQTPPI T33 w180mix 
    function [times, FID, FIDs] = TQTPPI_T33w180(obj)
            times = obj.deadtimeFID + obj.dwelltimeFID.*(0:obj.dataPoints-1);
            lengthAlpha = length(obj.alphas);
            tevo = obj.tevo0;
            tmixStart = 0.5 * 1e-3;
            
            nFIDs = obj.NumPhaseCycles * lengthAlpha;
    
            FIDs = zeros(nFIDs,obj.nSpins,length(times),2);
            for j = 1:obj.NumPhaseCycles
                for alphaIndex = 1:lengthAlpha
                    idx = (j-1)*lengthAlpha + alphaIndex ;
                    tmix = tmixStart +  (idx-1) * obj.tmixStep; 
                    tmixs(idx) = tmix;
                    idx = (j-1)*lengthAlpha + alphaIndex ;
                    alpha = obj.alphas(alphaIndex);
                    beta_p1 = alpha + 90;
                    beta_p2 = alpha - 90;
                    
                    FIDs_p1 = zeros(obj.nSpins,length(times));
                    FIDs_p2 = zeros(obj.nSpins,length(times));
                    
                    if obj.nSpins > 100
                        parfor idxSpin = 1:obj.nSpins 
                            w_Shift = obj.FreqShift + random(obj.wShift_dist);
                            
                            TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                            [~, ~, ~,~,Tmn_p1_Spin_relax] = PhaseCyclesDZ.general3Pulse_w180in23(TmnStart, alpha, beta_p1, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90, obj.flip180,0);
                            [~, ~, ~,~,Tmn_p2_Spin_relax] = PhaseCyclesDZ.general3Pulse_w180in23(TmnStart, alpha, beta_p2, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90, obj.flip180,0);
                 
                            FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                            FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                        end
                    else
                        for idxSpin = 1:obj.nSpins 
                            w_Shift = obj.FreqShift + random(obj.wShift_dist);
                            
                            TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                            [~, ~, ~,~,Tmn_p1_Spin_relax] = PhaseCyclesDZ.general3Pulse_w180in23(TmnStart, alpha, beta_p1, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90, obj.flip180,0);
                            [~, ~, ~,~,Tmn_p2_Spin_relax] = PhaseCyclesDZ.general3Pulse_w180in23(TmnStart, alpha, beta_p2, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90, obj.flip180,0);
                 
                            FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                            FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                        end
                    end
                    FIDs(idx,:,:,1) = FIDs_p1;
                    FIDs(idx,:,:,2) = FIDs_p2;
                    
                    fprintf('PhaceCycle %d/%d completed\n',idx,nFIDs);
                end
            end
                       
            FID_sumSpins = sum(FIDs,2);
            FID_sumFIDs = sum(FID_sumSpins,3);
            FID_sumBetas = FID_sumFIDs(:,:,:,1)+FID_sumFIDs(:,:,:,2);
    
            FID = getValues.switchReIm(FID_sumBetas);      
    
    
    
    
    end

    % TQTPPI T33 w180evo  
    function [tmixs, FID, FIDs] = evo180_TQT33no180(obj)
            times = obj.deadtimeFID + obj.dwelltimeFID.*(0:obj.dataPoints-1);
            lengthAlpha = length(obj.alphas);
            tevo = obj.tevo0;
            tmixStart = 0.150 * 1e-3;
            
            nFIDs = obj.NumPhaseCycles * lengthAlpha;
    
            FIDs = zeros(nFIDs,obj.nSpins,length(times),2);
            for j = 1:obj.NumPhaseCycles
                fprintf('Simulating Cycle No.: %f \n',j);
                for alphaIndex = 1:lengthAlpha
                    idx = (j-1)*lengthAlpha + alphaIndex ;
                    tmix = tmixStart +  (idx-1) * obj.tmixStep; 
                    tmixs(idx) = tmix;
                    idx = (j-1)*lengthAlpha + alphaIndex ;
                    alpha = obj.alphas(alphaIndex);
                    beta_p1 = alpha + 90;
                    beta_p2 = alpha - 90;
                    
                    FIDs_p1 = zeros(obj.nSpins,length(times));
                    FIDs_p2 = zeros(obj.nSpins,length(times));
                    
                    if obj.nSpins > 100

                        parfor idxSpin = 1:obj.nSpins 
                            w_Shift = obj.FreqShift + random(obj.wShift_dist);
                            
                            TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                            [~, ~, ~,~,Tmn_p1_Spin_relax] = PhaseCycle.general3Pulse_w180(TmnStart, alpha, beta_p1, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90, obj.flip180,0);
                            [~, ~, ~,~,Tmn_p2_Spin_relax] = PhaseCycle.general3Pulse_w180(TmnStart, alpha, beta_p2, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90, obj.flip180,0);
                 
                            FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                            FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                        end
                    else
                        for idxSpin = 1:obj.nSpins 
                            w_Shift = obj.FreqShift + random(obj.wShift_dist);
                            
                            TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                            [~, ~, ~,~,Tmn_p1_Spin_relax] = PhaseCycle.general3Pulse_w180(TmnStart, alpha, beta_p1, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90, obj.flip180,0);
                            [~, ~, ~,~,Tmn_p2_Spin_relax] = PhaseCycle.general3Pulse_w180(TmnStart, alpha, beta_p2, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90, obj.flip180,0);
                 
                            FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                            FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                        end
                    end
                    
                    FIDs(idx,:,:,1) = FIDs_p1;
                    FIDs(idx,:,:,2) = FIDs_p2;
                    
    %                 fprintf('PhaceCycle %d/%d completed\n',idx,nFIDs);
                end
            end
                       
            FID_sumSpins = sum(FIDs,2);
            FID_sumFIDs = sum(FID_sumSpins,3);
            FID_sumBetas = FID_sumFIDs(:,:,:,1)+FID_sumFIDs(:,:,:,2);
    
            FID = getValues.switchReIm(FID_sumBetas);      
    
    
    
    
    end

    % TQTPPI T33 no180 
    function [tmixs, FID, FIDs] = TQT33no180(obj)
            times = obj.deadtimeFID + obj.dwelltimeFID.*(0:obj.dataPoints-1);
            lengthAlpha = length(obj.alphas);
            tevo = obj.tevo0;
            tmixStart = 0.150 * 1e-3;
            
            nFIDs = obj.NumPhaseCycles * lengthAlpha;
    
            FIDs = zeros(nFIDs,obj.nSpins,length(times),2);
            for j = 1:obj.NumPhaseCycles
                fprintf('Simulating Cycle No.: %f \n',j);
                for alphaIndex = 1:lengthAlpha
                    idx = (j-1)*lengthAlpha + alphaIndex ;
                    tmix = tmixStart +  (idx-1) * obj.tmixStep; 
                    tmixs(idx) = tmix;
                    idx = (j-1)*lengthAlpha + alphaIndex ;
                    alpha = obj.alphas(alphaIndex);
                    beta_p1 = alpha + 90;
                    beta_p2 = alpha - 90;
                    
                    FIDs_p1 = zeros(obj.nSpins,length(times));
                    FIDs_p2 = zeros(obj.nSpins,length(times));
                    
                    if obj.nSpins > 100

                        parfor idxSpin = 1:obj.nSpins 
                            w_Shift = obj.FreqShift + random(obj.wShift_dist);
                            
                            TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                            [~, ~, ~,~,Tmn_p1_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha, beta_p1, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,0);
                            [~, ~, ~,~,Tmn_p2_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha, beta_p2, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,0);
                 
                            FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                            FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                        end
                    else
                        for idxSpin = 1:obj.nSpins 
                            w_Shift = obj.FreqShift + random(obj.wShift_dist);
                            
                            TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                            [~, ~, ~,~,Tmn_p1_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha, beta_p1, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90, 0);
                            [~, ~, ~,~,Tmn_p2_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha, beta_p2, tevo, tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90, 0);
                 
                            FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                            FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                        end
                    end
                    
                    FIDs(idx,:,:,1) = FIDs_p1;
                    FIDs(idx,:,:,2) = FIDs_p2;
                    
    %                 fprintf('PhaceCycle %d/%d completed\n',idx,nFIDs);
                end
            end
                       
            FID_sumSpins = sum(FIDs,2);
            FID_sumFIDs = sum(FID_sumSpins,3);
            FID_sumBetas = FID_sumFIDs(:,:,:,1)+FID_sumFIDs(:,:,:,2);
    
            FID = getValues.switchReIm(FID_sumBetas);      
    
    
    
    
    end

%% Other TQ sequences

    
    function [FIDs, FID, times] = sist2xNew_DZ(obj)
            times = obj.deadtimeFID +(0:obj.dataPoints-1) .* (obj.TR/obj.dataPoints);
            obj.alphas = [30:60:330];
            lengthAlpha = length(obj.alphas);
            tevo = obj.tevo0;
            nFIDs = obj.NumPhaseCycles * lengthAlpha;
            
            FIDs = zeros(nFIDs,obj.nSpins,length(times),2);
            for j = 1:obj.NumPhaseCycles
                for alphaIndex = 1:lengthAlpha
                    idx = (j-1)*lengthAlpha + alphaIndex ;
                    
                    alpha1 = obj.alphas(alphaIndex);
                    alpha2 = alpha1  + 180;
                    beta_p1 = alpha1 + 90;
                    beta_p2 = alpha2 - 90;
    
                    %phase_receiver = [0 180 0 180 0 180 ];
                    phase_receiver = [0 0 0 0 0 0];
                    FIDs_p1 = zeros(obj.nSpins,length(times));
                    FIDs_p2 = zeros(obj.nSpins,length(times));
                    %fprintf('AlphaIdx: %d \n', alphaIndex);
    
                    if obj.nSpins < 150
                        
                        for idxSpin = 1:obj.nSpins 
                            w_Shift = obj.FreqShift + random(obj.wShift_dist);
                            %flip90_shift = obj.flipShift + random(obj.flipShift_dist);
        %                     w_Shift = 0;
        
                            TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                            [~,~,~,~,Tmn_p1_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha1, beta_p1, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,phase_receiver(alphaIndex));
                            [~,~,~,~,Tmn_p2_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha2, beta_p2, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,phase_receiver(alphaIndex));
        
                            FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                            FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                        end
                    else
    
                        parfor idxSpin = 1:obj.nSpins 
                            w_Shift = obj.FreqShift + random(obj.wShift_dist);
        %                     w_Shift = 0;
        
                            TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                            [~,~,~,~,Tmn_p1_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha1, beta_p1, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,phase_receiver(alphaIndex));
                            [~,~,~,~,Tmn_p2_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha2, beta_p2, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,phase_receiver(alphaIndex));
        
                            FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                            FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                        end
                    end
                    FIDs(idx,:,:,1) = FIDs_p1;
                    FIDs(idx,:,:,2) = FIDs_p2;
                    
                    if mod(idx, 10) == 0
                        fprintf('PhaceCycle %d/%d completed\n',idx,nFIDs);
                    end
                end
            end
            FID_sumSpins = squeeze(sum(FIDs,2)); % spins are second axis
            FID_sumFIDs = mean(FID_sumSpins,1); % sum over the 6 cycles
            FID = FID_sumFIDs(:,:,1)+FID_sumFIDs(:,:,2); % add both
            %FID = getValues.switchReIm(FID);
        end



    function [FIDs, times] = Fleysher_DZ(obj)
            times = obj.deadtimeFID +(0:obj.dataPoints-1) .* (obj.TR/obj.dataPoints);
            obj.alphas = [0:60:300]; % 1st cycle
            lengthAlpha = length(obj.alphas);
            tevo = obj.tevo0;
            nFIDs = obj.NumPhaseCycles * lengthAlpha;
            
            FIDs = zeros(nFIDs,obj.nSpins,length(times),2);
            for j = 1:obj.NumPhaseCycles
                for alphaIndex = 1:lengthAlpha
                    idx = (j-1)*lengthAlpha + alphaIndex ;
                    
                    alpha1 = obj.alphas(alphaIndex);
                    alpha2 = alpha1  + 90; % 1st pulse, second cycle
                    beta_p1 = alpha1; % 2nd pulse, first cycle
                    beta_p2 = alpha2 - 90; % 2nd pulse, second cycle
    
                    %phase_receiver = [0 180 0 180 0 180 ];
                    phase_receiver = [0 0 0 0 0 0 ];
                    
                    FIDs_p1 = zeros(obj.nSpins,length(times));
                    FIDs_p2 = zeros(obj.nSpins,length(times));
                    %fprintf('AlphaIdx: %d \n', alphaIndex);
    
                    if obj.nSpins < 150
                        
                        for idxSpin = 1:obj.nSpins 
                            w_Shift = obj.FreqShift + random(obj.wShift_dist);
        %                     w_Shift = 0;
        
                            TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                            [~,~,~,~,Tmn_p1_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha1, beta_p1, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,phase_receiver(alphaIndex));
                            [~,~,~,~,Tmn_p2_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha2, beta_p2, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,phase_receiver(alphaIndex));
        
                            FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                            FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                        end
                    else
    
                        parfor idxSpin = 1:obj.nSpins 
                            w_Shift = obj.FreqShift + random(obj.wShift_dist);
        %                     w_Shift = 0;
        
                            TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                            [~,~,~,~,Tmn_p1_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha1, beta_p1, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,phase_receiver(alphaIndex));
                            [~,~,~,~,Tmn_p2_Spin_relax] = PhaseCycle.general3Pulse(TmnStart, alpha2, beta_p2, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,phase_receiver(alphaIndex));
        
                            FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                            FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
                        end
                    end
                    FIDs(idx,:,:,1) = FIDs_p1;
                    FIDs(idx,:,:,2) = FIDs_p2;
                    
                    if mod(idx, 10) == 0
                        fprintf('PhaceCycle %d/%d completed\n',idx,nFIDs);
                    end
                end
            end
            FID_sumSpins = squeeze(sum(FIDs,2)); % spins are second axis
            FID_sumFIDs = mean(FID_sumSpins,1); % sum over the 6 cycles
            FID = FID_sumFIDs(:,:,1)+FID_sumFIDs(:,:,2); % add both
            %FID = getValues.switchReIm(FID);
        end

function [FID_sumSpins, times] = boada_TQF_dz(obj)
    times = obj.deadtimeFID +(0:obj.dataPoints-1) .* (obj.TR/obj.dataPoints);
    lengthAlpha = length(obj.alphas);
    tevo = obj.tevo0;
    nFIDs = obj.NumPhaseCycles * lengthAlpha;

    Rx = [0 180 0 180 0 180];
    Rx2 = [180 0 180 0 180 0];

    
    FIDs = zeros(nFIDs,obj.nSpins,length(times),2);
    for j = 1:obj.NumPhaseCycles
        for alphaIndex = 1:lengthAlpha
            idx = (j-1)*lengthAlpha + alphaIndex ;
            alpha = obj.alphas(alphaIndex);
            beta_p1 = alpha + 90; % phase of second pulse first cycle
            alpha2 = alpha + 180;
            beta_p2 = beta_p1; % phase of second pulse second cycle
            
            FIDs_p1 = zeros(obj.nSpins,length(times));
            FIDs_p2 = zeros(obj.nSpins,length(times));
            for idxSpin = 1:obj.nSpins 
                w_Shift = obj.FreqShift + random(obj.wShift_dist);
%                     w_Shift = 0;
                
                TmnStart = TmnEvo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.Jen, w_Shift, obj.wShift_RMS, obj.wShift_FID);
                [~,~,~,~,Tmn_p1_Spin_relax] = PhaseCyclesDZ.general3Pulse(TmnStart, alpha, beta_p1, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,Rx(alphaIndex));
                [~,~,~,~,Tmn_p2_Spin_relax] = PhaseCyclesDZ.general3Pulse(TmnStart, alpha2, beta_p2, tevo, obj.tmix, obj.TR, obj.flip90, obj.flip90, obj.flip90,Rx2(alphaIndex));

                FIDs_p1(idxSpin,:) =  Tmn_p1_Spin_relax(1,5,1:end-1);
                FIDs_p2(idxSpin,:) =  Tmn_p2_Spin_relax(1,5,1:end-1); 
            end
            
            FIDs(idx,:,:,1) = FIDs_p1;
            FIDs(idx,:,:,2) = FIDs_p2;
%                 plot(-imag(sum(FIDs_p1,1)));
%                 fprintf('PhaceCycle %d/%d completed\n',idx,nFIDs);
        end
    end
    %first dimensional FID
%         FID_sumSpins = sum(FIDs,2);
%         FID_sumFIDs = sum(FID_sumSpins,1);
%         FID_sumFIDs = reshape(FID_sumFIDs,length(times),2);
    
    %whole FID for 1D spec
    FID_sumSpins = sum(FIDs,2);
    FID_sumFIDs = sum(FID_sumSpins,3);
%         FID_sumFIDs = reshape(FID_sumFIDs,nFIDs,2);
    
    %whole FID for 2D spec
%         FID_sumSpins = sum(FIDs,2);

    
    FIDp1 = getValues.switchReIm(FID_sumFIDs(:,:,:,1));
    FIDp2 = getValues.switchReIm(FID_sumFIDs(:,:,:,2));

      
end

end

end