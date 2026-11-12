function [A_Pseu, PHI2D_cell] = BEM2D(OFDMParams, OFDMConfig)

    N  = OFDMConfig.N;      % NUEVO: alias para abreviar
    CP = OFDMConfig.CP;     % NUEVO
    FB = OFDMParams.Frameblock;  % NUEVO

    %% Channel estimation 2D interpolation
    %DPSS BEM configuration to time - delay domain with RU 
    seq_lengthd = OFDMConfig.CP+1;%50;       %Infered from power delay profile, can be extended to CP
    %Normalized time-bandwidth utilization
    Xlength=OFDMParams.Nu+1;
    time_halfbandwidth = ((Xlength)/(2*(OFDMConfig.N)))*seq_lengthd;
    %time_halfbandwidth = (OFDMconfig.cP-1)/2; %Time halfbandwidth +1
    num_seq=ceil(2*time_halfbandwidth)+2; %Number of functions tenia 1,puse 3 y deje 2
    [ExpFuncdelay,~] = dpss(seq_lengthd,time_halfbandwidth,num_seq);
    
    %Time basis for complete subframe duration, OFDM sub block-resolution 
    seq_length = (OFDMConfig.N+OFDMConfig.CP)*(OFDMParams.Frameblock);
    time_halfbandwidthtime = (2*((OFDMConfig.fd)*(OFDMConfig.N+OFDMConfig.CP)*...
       (OFDMParams.Frameblock)*OFDMConfig.Ts))/2;
    NumSeqTime=ceil(2*time_halfbandwidthtime)+3;%2 3
    [TimeSequences,~] = dpss(seq_length,time_halfbandwidthtime,NumSeqTime);
    
    %Initialize and fill 2D BEMs
    TwoDbasisNum=num_seq*NumSeqTime;

    % ---- ORIGINAL: array 3D completo en memoria (~545 MB) ----
    % Phi2D=zeros(OFDMConfig.CP+1,(OFDMConfig.N+OFDMConfig.CP)*(OFDMParams.Frameblock),TwoDbasisNum);
    % for j=1:NumSeqTime
    %     for i=1:num_seq
    %         Phi2D(:,:,(j-1)*num_seq+i)=ExpFuncdelay(:,i)*transpose(TimeSequences(:,j));
    %     end
    % end
    % ---- NUEVO: no se construye Phi2D 3D; cada base se calcula en su iteración ----
    block_len = seq_length / FB;   % = N+CP, longitud de cada bloque

    %% Phi2D without CP
    % ---- ORIGINAL: PHI2D_cell se llenaba a partir de Phi2D ya precalculado ----
    % iter2 = 0;
    % PHI2D_cell = cell(TwoDbasisNum,(OFDMParams.Frameblock));
    % H = zeros(OFDMConfig.N,OFDMConfig.N);
    % for i=1:TwoDbasisNum
    %     temp = mat2cell(Phi2D(:,:,i),numel(Phi2D(:,1,1)),repmat(numel(Phi2D(1,:,1))/(OFDMParams.Frameblock),1,(OFDMParams.Frameblock)));
    %     for ii=1:length(temp)
    %         v = temp{1,ii};
    %         v = v(:,OFDMConfig.CP+1:end);
    %         for iii = 1:OFDMConfig.N
    %             v_n = (v(:,iii)).';
    %             v_n_ZP = [v_n zeros(1,OFDMConfig.N-numel(v_n))];
    %             ZP2 = [v_n_ZP(1) fliplr(v_n_ZP(2:end))];
    %             ZP_1 = circshift(ZP2,iter2);
    %             H(iii,:) = ZP_1;
    %             iter2 = iter2 + 1;
    %         end
    %         iter2 = 0;
    %         H_Toep_Freq = ifft(fft(H,[],1),[],2);
    %         H_Toep_Freq = fftshift(H_Toep_Freq,1);
    %         H_Toep_Freq = ifftshift(H_Toep_Freq,2);
    %         PHI2D_cell{i,ii} = H_Toep_Freq;
    %     end
    %     i
    % end

    % ---- NUEVO: se llena PHI2D_cell base por base, sin Phi2D intermedio.
    %      Adicionalmente se calcula RowSums mientras se recorre, para no
    %      tener que releer PHI2D_cell más adelante al construir A_Pseu.
    PHI2D_cell = cell(TwoDbasisNum, FB);
    RowSums    = zeros(N, FB, TwoDbasisNum);   % ~6.5 MB, reemplaza al bloque PHI muerto
    H          = zeros(N, N);
    onevec     = ones(N,1);

    for k = 1:TwoDbasisNum
        % índice (i,j) equivalente al bucle doble original j=1:NumSeqTime, i=1:num_seq
        j_time  = ceil(k / num_seq);
        i_delay = k - (j_time-1)*num_seq;

        % Base 2D para esta k: (CP+1) x seq_length   [equivalente a Phi2D(:,:,k)]
        Phi2D_k = ExpFuncdelay(:,i_delay) * TimeSequences(:,j_time).';

        for ii = 1:FB
            % reemplaza mat2cell + temp{1,ii}
            v = Phi2D_k(:, (ii-1)*block_len+1 : ii*block_len);
            v = v(:, CP+1:end);                 % quitar CP

            for iii = 1:N
                v_n = v(:,iii).';
                v_n_ZP = [v_n zeros(1, N-numel(v_n))];
                ZP2  = [v_n_ZP(1) fliplr(v_n_ZP(2:end))];
                % iter2 en el original equivale siempre a (iii-1) dentro del bucle
                H(iii,:) = circshift(ZP2, iii-1);
            end

            H_Toep_Freq = ifft(fft(H,[],1),[],2);
            H_Toep_Freq = fftshift(H_Toep_Freq,1);
            H_Toep_Freq = ifftshift(H_Toep_Freq,2);

            PHI2D_cell{k,ii} = H_Toep_Freq;
            RowSums(:,ii,k)  = H_Toep_Freq * onevec;   % NUEVO: guarda la suma que A necesita
        end
        clear Phi2D_k
    end
    clear H

    % ---- ORIGINAL: bloque PHI (dead code, nunca se usaba después) ----
    % suma = 0;
    % PHI = cell(1,OFDMParams.Frameblock);
    % for dann = 1:OFDMParams.Frameblock
    %     for dan = 1:TwoDbasisNum 
    %         suma = suma + PHI2D_cell{dan,dann};
    %     end
    %     PHI{1,dann} = suma; 
    %     suma = 0;
    % end

    % ---- ORIGINAL: pipeline PHI2D_Vec -> cell2mat -> PHI2D_COMP (dead code) ----
    % PHI2D_cell2 = transpose(PHI2D_cell);
    % PHI2D_Vec = cell(OFDMParams.Frameblock,TwoDbasisNum);
    % for a = 1:TwoDbasisNum
    %     for aa = 1:OFDMParams.Frameblock
    %         auxiliar2 = PHI2D_cell2{aa,a};
    %         PHI2D_Vec{aa,a} = auxiliar2(:);
    %     end
    % end 
    % PHI2D_Vec = (cell2mat(PHI2D_Vec));
    % PHI2D_COMP = mat2cell(PHI2D_Vec,repmat(OFDMConfig.N^2,size(PHI2D_Vec,1)/OFDMConfig.N^2,1),TwoDbasisNum);

    %Receiver model to compute the weigths and estimation matrix for one and
    %other symbols
    % ---- ORIGINAL: cell array A + cell2mat ----
    % A = cell(TwoDbasisNum,OFDMParams.Frameblock);
    % for j=1:OFDMParams.Frameblock
    %     for jj=1:TwoDbasisNum 
    %         VecPhi2D = PHI2D_cell2{j,jj};
    %         A{jj,j}=VecPhi2D*ones(OFDMConfig.N,1);
    %     end  
    % end
    % A = transpose(A);
    % A2 = cell2mat(A);
    % A2 = A2(OFDMParams.RUPilotFBIndexesVec,:);
    % A_INV = (A2'*A2);
    % clear A;
    % A_Pseu = A_INV\A2';
    % clear A2;

    % ---- NUEVO: A2 se construye directamente desde RowSums ----
    %      Layout idéntico al original: (FB*N) x TwoDbasisNum, bloques por frame.
    A2 = zeros(FB*N, TwoDbasisNum);
    for j = 1:FB
        for jj = 1:TwoDbasisNum
            A2((j-1)*N+1 : j*N, jj) = RowSums(:, j, jj);
        end
    end
    clear RowSums

    A2 = A2(OFDMParams.RUPilotFBIndexesVec, :);
    A_INV  = A2' * A2;
    A_Pseu = A_INV \ A2';
    clear A2 A_INV

end