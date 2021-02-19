% Script made for performing general tests.

function runall_segmentation(path,pathOut)

    folders = dir(path)
   
    for sdr = 3:numel(folders)
        
        subfolder = folders(sdr).name;
        
        pngPath = sprintf('%s%s/*.png',path,subfolder);
        
        lista = dir(pngPath);
        
        tamLista = size(lista);
        
        for i =1:tamLista
        
            segPath = sprintf('%s%s/*.png',pathOut,subfolder);
            
            imgPath = sprintf('%s%s/%s',path,subfolder,lista(i).name);

            segPath = sprintf('%s%s/%s_bin.png',pathOut,subfolder,lista(i).name(1:end-4));
            
            fprintf('%d - %s - %s\n',i,imgPath,segPath);
            
            imgRaw = double(imread(imgPath));
            [height,width,channels] = size(imgRaw);


            maior = max(max(imgRaw(:,:,2)));
            if(maior > 255)
                imgRaw = uint8(imgRaw / 256);
            end

            % Take eigenvalues
            rawHSV = rgb2hsv(uint8(imgRaw));
            satur = rawHSV(:,:,2);
            w = fspecial('gaussian',15,2);
            gausssatur = imfilter(satur,w);


            binMask = gausssatur > mean(mean(gausssatur))+0.05;
            algaeSamplingMask = imopen(binMask,strel('disk',2));
            algaeSamplingMask = imerode(algaeSamplingMask,strel('disk',2));

            [height,width,canais] = size(imgRaw);

            n=1;
            nn = 1;
            for k = 1:height
                for j = 1:width

                    if algaeSamplingMask(k,j) == 1

                        Reg1(n,1) = double(imgRaw(k,j,1));
                        Reg1(n,2) = double(imgRaw(k,j,2));
                        Reg1(n,3) = double(imgRaw(k,j,3));
                        n = n+1;

                    else

                        largeReg2(nn,1) = double(imgRaw(k,j,1));
                        largeReg2(nn,2) = double(imgRaw(k,j,2));
                        largeReg2(nn,3) = double(imgRaw(k,j,3));
                        nn = nn+1;
                    end

                end
            end

            factor = max(round(nn/20),round(nn/50));
            Reg2(1:factor,1:3) = 0;
            Reg2(:,1) = largeReg2(1:factor,1);
            Reg2(:,2) = largeReg2(1:factor,2);
            Reg2(:,3) = largeReg2(1:factor,3);

            mask = imgRaw;

            nnn = 1;
            for k = 1:height
                for j = 1:width

                    if binMask(k,j) == 1
                        mask(k,j,1) = 255;
                        mask(k,j,2) = 0;
                        mask(k,j,3) = 0;
                    else
                        if nnn < factor
                            mask(k,j,1) = 0;
                            mask(k,j,2) = 255;
                            mask(k,j,3) = 0;
                        end
                    end
                    nnn = nnn+1;
                end
            end
            
            u = levelset(imgPath,0.1,[],0.05,0.0001,1,0,1800,segPath,Reg1,Reg2,0,segPath,lista(i).name(1:end-4),1);

            close all;
        
            clear Reg1;
            clear Reg2;
            clear largeReg2;
            clear largeReg1;
        
        end
    end
end