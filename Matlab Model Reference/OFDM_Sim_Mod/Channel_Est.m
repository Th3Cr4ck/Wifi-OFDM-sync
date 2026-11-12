function G = Channel_Est(OFDMParams,OFDMConfig,A_Pseu,PHI2D_cell,RxFrecFrame)
    RxFrecFrame_Vec = RxFrecFrame(:);
    count = 1;
    has_pilots = zeros(OFDMParams.M,1);
    for i = 1:OFDMParams.MRU
        for ii = 1:OFDMParams.Frameblock
            if (OFDMParams.PerSymbol{ii, 1}.hasPilots == 1)
                has_pilots(count,1) = 1;
                count = count + 1;
            else
                count = count + 1;
            end
        end
    end
    
    numb_pilot_symb = numel(find(has_pilots == 1));
    x_p_vec = OFDMParams.TxPilots;
    y_p_vec = RxFrecFrame_Vec(OFDMParams.PilotIndexesVec);
    i_pilots = y_p_vec./x_p_vec;
  
    pilots_val_y_p = reshape(i_pilots,[numel(OFDMParams.PerSymbol{1,1}.Pilot),numb_pilot_symb]);
    
    %% Wheight Estimation 
    switch OFDMConfig.MRU
        case 1
            pilots_val_y_p1 = pilots_val_y_p(:,1:4);

            Rho_Weigths = A_Pseu*pilots_val_y_p1(:); 
        case 2
            pilots_val_y_p1 = pilots_val_y_p(:,1:4);
            pilots_val_y_p2 = pilots_val_y_p(:,5:8);
            
            Rho_Weigths = A_Pseu*pilots_val_y_p1(:); 
            Rho_Weigths2 = A_Pseu*pilots_val_y_p2(:);
        case 3
            pilots_val_y_p1 = pilots_val_y_p(:,1:4);
            pilots_val_y_p2 = pilots_val_y_p(:,5:8);
            pilots_val_y_p3 = pilots_val_y_p(:,9:12);

            Rho_Weigths = A_Pseu*pilots_val_y_p1(:); 
            Rho_Weigths2 = A_Pseu*pilots_val_y_p2(:);
            Rho_Weigths3 = A_Pseu*pilots_val_y_p3(:); 
        case 4
            pilots_val_y_p1 = pilots_val_y_p(:,1:4);
            pilots_val_y_p2 = pilots_val_y_p(:,5:8);
            pilots_val_y_p3 = pilots_val_y_p(:,9:12);
            pilots_val_y_p4 = pilots_val_y_p(:,13:16);

            Rho_Weigths = A_Pseu*pilots_val_y_p1(:); 
            Rho_Weigths2 = A_Pseu*pilots_val_y_p2(:);
            Rho_Weigths3 = A_Pseu*pilots_val_y_p3(:); 
            Rho_Weigths4 = A_Pseu*pilots_val_y_p4(:);
        case 5
            pilots_val_y_p1 = pilots_val_y_p(:,1:4);
            pilots_val_y_p2 = pilots_val_y_p(:,5:8);
            pilots_val_y_p3 = pilots_val_y_p(:,9:12);
            pilots_val_y_p4 = pilots_val_y_p(:,13:16);
            pilots_val_y_p5 = pilots_val_y_p(:,17:20);
            
            Rho_Weigths = A_Pseu*pilots_val_y_p1(:); 
            Rho_Weigths2 = A_Pseu*pilots_val_y_p2(:);
            Rho_Weigths3 = A_Pseu*pilots_val_y_p3(:); 
            Rho_Weigths4 = A_Pseu*pilots_val_y_p4(:);
            Rho_Weigths5 = A_Pseu*pilots_val_y_p5(:); 
        case 6
            pilots_val_y_p1 = pilots_val_y_p(:,1:4);
            pilots_val_y_p2 = pilots_val_y_p(:,5:8);
            pilots_val_y_p3 = pilots_val_y_p(:,9:12);
            pilots_val_y_p4 = pilots_val_y_p(:,13:16);
            pilots_val_y_p5 = pilots_val_y_p(:,17:20);
            pilots_val_y_p6 = pilots_val_y_p(:,21:24);

            Rho_Weigths = A_Pseu*pilots_val_y_p1(:); 
            Rho_Weigths2 = A_Pseu*pilots_val_y_p2(:);
            Rho_Weigths3 = A_Pseu*pilots_val_y_p3(:); 
            Rho_Weigths4 = A_Pseu*pilots_val_y_p4(:);
            Rho_Weigths5 = A_Pseu*pilots_val_y_p5(:); 
            Rho_Weigths6 = A_Pseu*pilots_val_y_p6(:);
        case 7
            pilots_val_y_p1 = pilots_val_y_p(:,1:4);
            pilots_val_y_p2 = pilots_val_y_p(:,5:8);
            pilots_val_y_p3 = pilots_val_y_p(:,9:12);
            pilots_val_y_p4 = pilots_val_y_p(:,13:16);
            pilots_val_y_p5 = pilots_val_y_p(:,17:20);
            pilots_val_y_p6 = pilots_val_y_p(:,21:24);
            pilots_val_y_p7 = pilots_val_y_p(:,25:28);

            Rho_Weigths = A_Pseu*pilots_val_y_p1(:); 
            Rho_Weigths2 = A_Pseu*pilots_val_y_p2(:);
            Rho_Weigths3 = A_Pseu*pilots_val_y_p3(:); 
            Rho_Weigths4 = A_Pseu*pilots_val_y_p4(:);
            Rho_Weigths5 = A_Pseu*pilots_val_y_p5(:); 
            Rho_Weigths6 = A_Pseu*pilots_val_y_p6(:);
            Rho_Weigths7 = A_Pseu*pilots_val_y_p7(:); 
        case 8
             pilots_val_y_p1 = pilots_val_y_p(:,1:4);
            pilots_val_y_p2 = pilots_val_y_p(:,5:8);
            pilots_val_y_p3 = pilots_val_y_p(:,9:12);
            pilots_val_y_p4 = pilots_val_y_p(:,13:16);
            pilots_val_y_p5 = pilots_val_y_p(:,17:20);
            pilots_val_y_p6 = pilots_val_y_p(:,21:24);
            pilots_val_y_p7 = pilots_val_y_p(:,25:28);
            pilots_val_y_p8 = pilots_val_y_p(:,29:32);

            Rho_Weigths = A_Pseu*pilots_val_y_p1(:); 
            Rho_Weigths2 = A_Pseu*pilots_val_y_p2(:);
            Rho_Weigths3 = A_Pseu*pilots_val_y_p3(:); 
            Rho_Weigths4 = A_Pseu*pilots_val_y_p4(:);
            Rho_Weigths5 = A_Pseu*pilots_val_y_p5(:); 
            Rho_Weigths6 = A_Pseu*pilots_val_y_p6(:);
            Rho_Weigths7 = A_Pseu*pilots_val_y_p7(:); 
            Rho_Weigths8 = A_Pseu*pilots_val_y_p8(:);
    end

    %% Channel Estimation 
    PHI2D_cell = transpose(PHI2D_cell);
    I = numel(Rho_Weigths);
    
    switch OFDMConfig.MRU
        case 1
            acumm = 0;
            ChannelCoeffs1 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths(ii,1);
                end
    
            %ChannelCoeffs1{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs1{i,1} = acumm;
            acumm = 0;
            end

            G = ChannelCoeffs1;
            clear ChannelCoeffs1; clear Rho_Weigths; clear PHI2D_cell;
        case 2
            acumm = 0;
            ChannelCoeffs1 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths(ii,1);
                end
            
            %ChannelCoeffs1{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs1{i,1} = acumm;
            acumm = 0;
            end

            acumm = 0;
            ChannelCoeffs2 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths2(ii,1);
                end
            
               %ChannelCoeffs2{i,1} = fftshift(fftshift(acumm,1),2);
               ChannelCoeffs2{i,1} = acumm;
               acumm = 0;
            end
            
            G = cat(1,ChannelCoeffs1,ChannelCoeffs2);
            clear ChannelCoeffs1; clear Rho_Weigths; clear PHI2D_cell;
            clear ChannelCoeffs2; clear Rho_Weigths2;
        
        case 3
            acumm = 0;
            ChannelCoeffs1 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths(ii,1);
                end
            
            %ChannelCoeffs1{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs1{i,1} = acumm;
            acumm = 0;
            end

            acumm = 0;
            ChannelCoeffs2 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths2(ii,1);
                end
        
               %ChannelCoeffs2{i,1} = fftshift(fftshift(acumm,1),2);
               ChannelCoeffs2{i,1} = acumm;
               acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs3 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths3(ii,1);
                end
        
            %ChannelCoeffs3{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs3{i,1} = acumm;
            acumm = 0;
            end
            G = cat(1,ChannelCoeffs1,ChannelCoeffs2,ChannelCoeffs3);
            clear ChannelCoeffs1; clear Rho_Weigths; clear PHI2D_cell;
            clear ChannelCoeffs2; clear Rho_Weigths2;
            clear ChannelCoeffs3; clear Rho_Weigths3;

        case 4
            acumm = 0;
            ChannelCoeffs1 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths(ii,1);
                end
            
            %ChannelCoeffs1{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs1{i,1} = acumm;
            acumm = 0;
            end

            acumm = 0;
            ChannelCoeffs2 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths2(ii,1);
                end
        
               %ChannelCoeffs2{i,1} = fftshift(fftshift(acumm,1),2);
               ChannelCoeffs2{i,1} = acumm;
               acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs3 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths3(ii,1);
                end
        
            %ChannelCoeffs3{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs3{i,1} = acumm;
            acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs4 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths4(ii,1);
                end
        
            %ChannelCoeffs4{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs4{i,1} = acumm;
            acumm = 0;
            end
            G = cat(1,ChannelCoeffs1,ChannelCoeffs2,ChannelCoeffs3,ChannelCoeffs4);
            clear ChannelCoeffs1; clear Rho_Weigths; clear PHI2D_cell;
            clear ChannelCoeffs2; clear Rho_Weigths2;
            clear ChannelCoeffs3; clear Rho_Weigths3;
            clear ChannelCoeffs4; clear Rho_Weigths4;

        case 5
            acumm = 0;
            ChannelCoeffs1 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths(ii,1);
                end
            
            %ChannelCoeffs1{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs1{i,1} = acumm;
            acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs2 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths2(ii,1);
                end

               %ChannelCoeffs2{i,1} = fftshift(fftshift(acumm,1),2);
               ChannelCoeffs2{i,1} = acumm;
               acumm = 0;
            end

            acumm = 0;
            ChannelCoeffs3 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths3(ii,1);
                end

            %ChannelCoeffs3{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs3{i,1} = acumm;
            acumm = 0;
            end

            acumm = 0;
            ChannelCoeffs4 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths4(ii,1);
                end

            %ChannelCoeffs4{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs4{i,1} = acumm;
            acumm = 0;
            end

            acumm = 0;
            ChannelCoeffs5 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths5(ii,1);
                end

            %ChannelCoeffs5{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs5{i,1} = acumm;
            acumm = 0;
            end
            G = cat(1,ChannelCoeffs1,ChannelCoeffs2,ChannelCoeffs3,ChannelCoeffs4,ChannelCoeffs5);
            clear ChannelCoeffs1; clear Rho_Weigths; clear PHI2D_cell;
            clear ChannelCoeffs2; clear Rho_Weigths2;
            clear ChannelCoeffs3; clear Rho_Weigths3;
            clear ChannelCoeffs4; clear Rho_Weigths4;
            clear ChannelCoeffs5; clear Rho_Weigths5;
        case 6
             acumm = 0;
            ChannelCoeffs1 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths(ii,1);
                end
            
            %ChannelCoeffs1{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs1{i,1} = acumm;
            acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs2 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths2(ii,1);
                end

               %ChannelCoeffs2{i,1} = fftshift(fftshift(acumm,1),2);
               ChannelCoeffs2{i,1} = acumm;
               acumm = 0;
            end

            acumm = 0;
            ChannelCoeffs3 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths3(ii,1);
                end

            %ChannelCoeffs3{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs3{i,1} = acumm;
            acumm = 0;
            end

            acumm = 0;
            ChannelCoeffs4 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths4(ii,1);
                end

            %ChannelCoeffs4{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs4{i,1} = acumm;
            acumm = 0;
            end

            acumm = 0;
            ChannelCoeffs5 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths5(ii,1);
                end

            %ChannelCoeffs5{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs5{i,1} = acumm;
            acumm = 0;
            end

            acumm = 0;
            ChannelCoeffs6 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths6(ii,1);
                end

            %ChannelCoeffs6{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs6{i,1} = acumm;
            acumm = 0;
            end
            G = cat(1,ChannelCoeffs1,ChannelCoeffs2,ChannelCoeffs3,ChannelCoeffs4,ChannelCoeffs5, ...
            ChannelCoeffs6);
            clear ChannelCoeffs1; clear Rho_Weigths; clear PHI2D_cell;
            clear ChannelCoeffs2; clear Rho_Weigths2;
            clear ChannelCoeffs3; clear Rho_Weigths3;
            clear ChannelCoeffs4; clear Rho_Weigths4;
            clear ChannelCoeffs5; clear Rho_Weigths5;
            clear ChannelCoeffs6; clear Rho_Weigths6;
        case 7 
            acumm = 0;
            ChannelCoeffs1 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths(ii,1);
                end
            
            %ChannelCoeffs1{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs1{i,1} = acumm;
            acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs2 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths2(ii,1);
                end
        
               %ChannelCoeffs2{i,1} = fftshift(fftshift(acumm,1),2);
               ChannelCoeffs2{i,1} = acumm;
               acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs3 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths3(ii,1);
                end
        
            %ChannelCoeffs3{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs3{i,1} = acumm;
            acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs4 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths4(ii,1);
                end
        
            %ChannelCoeffs4{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs4{i,1} = acumm;
            acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs5 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths5(ii,1);
                end
        
            %ChannelCoeffs5{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs5{i,1} = acumm;
            acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs6 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths6(ii,1);
                end
        
            %ChannelCoeffs6{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs6{i,1} = acumm;
            acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs7 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths7(ii,1);
                end
        
            %ChannelCoeffs7{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs7{i,1} = acumm;
            acumm = 0;
            end
            G = cat(1,ChannelCoeffs1,ChannelCoeffs2,ChannelCoeffs3,ChannelCoeffs4,ChannelCoeffs5, ...
            ChannelCoeffs6,ChannelCoeffs7);
            clear ChannelCoeffs1; clear Rho_Weigths; clear PHI2D_cell;
            clear ChannelCoeffs2; clear Rho_Weigths2;
            clear ChannelCoeffs3; clear Rho_Weigths3;
            clear ChannelCoeffs4; clear Rho_Weigths4;
            clear ChannelCoeffs5; clear Rho_Weigths5;
            clear ChannelCoeffs6; clear Rho_Weigths6;
            clear ChannelCoeffs7; clear Rho_Weigths7;
        case 8
            acumm = 0;
            ChannelCoeffs1 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths(ii,1);
                end
            
            %ChannelCoeffs1{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs1{i,1} = acumm;
            acumm = 0;
            end
        
            acumm = 0;
            ChannelCoeffs2 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths2(ii,1);
                end
    
               %ChannelCoeffs2{i,1} = fftshift(fftshift(acumm,1),2);
               ChannelCoeffs2{i,1} = acumm;
               acumm = 0;
            end
    
            acumm = 0;
            ChannelCoeffs3 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths3(ii,1);
                end
    
            %ChannelCoeffs3{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs3{i,1} = acumm;
            acumm = 0;
            end
    
            acumm = 0;
            ChannelCoeffs4 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths4(ii,1);
                end
    
            %ChannelCoeffs4{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs4{i,1} = acumm;
            acumm = 0;
            end
    
            acumm = 0;
            ChannelCoeffs5 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths5(ii,1);
                end
    
            %ChannelCoeffs5{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs5{i,1} = acumm;
            acumm = 0;
            end
    
            acumm = 0;
            ChannelCoeffs6 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths6(ii,1);
                end
    
            %ChannelCoeffs6{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs6{i,1} = acumm;
            acumm = 0;
            end
    
            acumm = 0;
            ChannelCoeffs7 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths7(ii,1);
                end
    
            %ChannelCoeffs7{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs7{i,1} = acumm;
            acumm = 0;
            end
    
            acumm = 0; 
            ChannelCoeffs8 = cell(OFDMParams.Frameblock,1);
            for i = 1:OFDMParams.Frameblock
                for ii = 1:I
                    acumm  = acumm + PHI2D_cell{i,ii}*Rho_Weigths8(ii,1);
                end
    
            %ChannelCoeffs8{i,1} = fftshift(fftshift(acumm,1),2);
            ChannelCoeffs8{i,1} = acumm;
            acumm = 0;
            end
            G = cat(1,ChannelCoeffs1,ChannelCoeffs2,ChannelCoeffs3,ChannelCoeffs4,ChannelCoeffs5, ...
            ChannelCoeffs6,ChannelCoeffs7,ChannelCoeffs8);
            clear ChannelCoeffs1; clear Rho_Weigths; clear PHI2D_cell;
            clear ChannelCoeffs2; clear Rho_Weigths2;
            clear ChannelCoeffs3; clear Rho_Weigths3;
            clear ChannelCoeffs4; clear Rho_Weigths4;
            clear ChannelCoeffs5; clear Rho_Weigths5;
            clear ChannelCoeffs6; clear Rho_Weigths6;
            clear ChannelCoeffs7; clear Rho_Weigths7;
            clear ChannelCoeffs8; clear Rho_Weigths8;
    end
end 