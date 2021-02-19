% u = regionCompetition('087_f2_ex3_rszd.png',0.1,[98;108;71;86;50;83;23;113],0.5,0.0001,1,0.0000005,2000,'t05');

function u = levelset(path,dt,amostra,LAMBDA,VU,EPSILON,BETA,T,outputPath,amost1,amost2,semiSup,segPath,imgName,modelOption)

    %PARAMETROS
   
    %Le imagem do disco
    imgRaw = imread(path);
    img = imresize(imgRaw,0.5);
    [altura,largura,canais] = size(img);
    
    I = padarray((double(img)), [1 1], 'replicate');
    precomputed =0;
    if precomputed == 0
        
        if(semiSup == 0)
            Reg1 = amost1;
            Reg2 = amost2;
        else
            Reg1 = double(img(amostra(1):amostra(2),amostra(3):amostra(4),:));
            Reg2 = double(img(amostra(5):amostra(6),amostra(7):amostra(8),:));
        end

        M1(1) = mean(Reg1(:,1));
        M1(2) = mean(Reg1(:,2));
        M1(3) = mean(Reg1(:,3));

        M2(1) = mean(Reg2(:,1));
        M2(2) = mean(Reg2(:,2));
        M2(3) = mean(Reg2(:,3));

        soma = 0;

        for i = 1:altura+2
            for j = 1:largura+2
                x(1) = I(i,j,1);
                x(2) = I(i,j,2);
                x(3) = I(i,j,3);
                soma = soma + ( (  (double(x) - M1 )  )' * ( double(x) - M1 ) );
            end
        end

        Cov1(:,:) = (1/(altura*largura-1))*soma;

        soma = 0;

        for i = 1:altura+2
            for j = 1:largura+2
                x(1) = I(i,j,1);
                x(2) = I(i,j,2);
                x(3) = I(i,j,3);
                soma = soma + ( (  (double(x) - M2 )  )' * ( double(x) - M2 ) );
            end
        end

        Cov2(:,:) = (1/(altura*largura-1))*soma;

    %         x = [-2000:.5:2000];
    %         norm = normpdf(x,M1(2),Cov1(2,2));
    %         norm2 = normpdf(x,M2(2),Cov2(2,2));
    %         figure,plot(x,norm,'r',x,norm2,'g','LineWidth',3);
    %         str = sprintf('Gaussian distributions for the green channel: Red patches: algae regions; Green patch: background');
    %         filename = sprintf('gaussLSM_%s.png','087_6dia_2.png');
    %         title(str,'FontSize',12,'FontWeight','bold');
    %         set(gca,'FontSize',12,'FontWeight','bold');
    %         grid on;
    %         print(filename,'-dpng');
    %         pause;

        for i = 1:altura+2
            for j = 1:largura+2
                x(1) = I(i,j,1);
                x(2) = I(i,j,2);
                x(3) = I(i,j,3);
                %P1(i,j) = multiGaussian(x,M1,Cov1);
                P1(i,j)=comp_gauss_dens_val(M1',Cov1,x');

                %P2(i,j) = multiGaussian(x,M2,Cov2);
                P2(i,j)=comp_gauss_dens_val(M2',Cov2,x');
            end
        end

        pathcomp = sprintf('competition/%s.mat',imgName);
        competition = log10((P2+ 1e-100)./ (P1 + 1e-100));
        save(pathcomp,'competition','M1','M2','Cov1','Cov2','P1','P2');
    else
        pathcomp = sprintf('competition/%s.mat',imgName);
        load(pathcomp,'competition','M1','M2','P1','P2');
    end
    
    I = rgb2gray(I);
    

    u = zeros(altura+2,largura+2);
    u = double(u);
    IGrad = u;
    delta = u;
    un = u;
    
    c1 = floor(altura/2);
    c2 = floor(largura/2);
    
    u(:,:) = +4;
    u(25:altura-20,25:largura-20) = -4;
    u(25,25:largura-20) = 0;
    u(altura-20,25:largura-20) = 0;
    u(25:altura-20,25) = 0;
    u(25:altura-20,largura-20) = 0;

    for i = 2:altura+1
        for j = 2:largura+1
             dudx = (I(i+1,j) - I(i-1,j))/2;
             dudy = (I(i,j+1) - I(i,j-1))/2;
             IGrad(i,j) = dudx^2 + dudy^2;
        end
    end
    
     g = 1 ./ (1 + BETA*IGrad);
    %g(1:altura+2,1:largura+2)=1;
    for i = 2:altura+1
        for j= 2:largura+1
            gx(i,j) = (g(i+1,j) - g(i-1,j))/2;
            gy(i,j) = (g(i,j+1) - g(i,j-1))/2;
        end
    end
    
    curva = bwboundaries(u);
    
    for i = 2:altura+1
        for j = 2:largura+1
            d = distCurva(i,j,curva{1,1});
            if d == 0
                u(i,j) = 0;
            else
                if u(i,j) == 4
                    u(i,j) = d;
                else
                    u(i,j) = -d;
                end
            end
        end
    end

    
    
    for t = 1 : T
        %Atualiza competition
        if modelOption == 2
            m1 = sum(sum(u(2:altura+1,2:largura+1)));
            m2 = sum(sum(1-u(2:altura+1,2:largura+1)));
            
            for c =1:3
                M1(1) = sum(sum(I(2:altura+1,2:largura+1,1).*u(2:altura+1,2:largura+1)))/m1;
                M1(2) = sum(sum(I(2:altura+1,2:largura+1,2).*u(2:altura+1,2:largura+1)))/m1;
                M1(3) = sum(sum(I(2:altura+1,2:largura+1,3).*u(2:altura+1,2:largura+1)))/m1;
                
                M2(1) = sum(sum(I(2:altura+1,2:largura+1,1).*(1-u(2:altura+1,2:largura+1))))/m2;
                M2(2) = sum(sum(I(2:altura+1,2:largura+1,2).*(1-u(2:altura+1,2:largura+1))))/m2;
                M2(3) = sum(sum(I(2:altura+1,2:largura+1,3).*(1-u(2:altura+1,2:largura+1))))/m2;
            end
            
            
            
            %competition2 = (((I-M2).^2)) - (((I-M1).^2));
        end
    
        %Atualiza evolucao das curvas
        for i = 2:altura+1
            for j= 2:largura+1
                
                
                %Calcula derivadas aproximadas de u
                ux = (u(i+1,j) - u(i-1,j))/2;
                uy = (u(i,j+1) - u(i,j-1))/2;
                uxx = (u(i+1,j) - 2*u(i,j) + u(i-1,j));
                uyy = (u(i,j+1) - 2*u(i,j) + u(i,j-1));
                uxy = (u(i+1,j+1) - u(i+1,j-1) - u(i-1,j+1) + u(i-1,j-1))/4;
                
                u2x = ux*ux;
                u2y = uy*uy;
                mag = sqrt(u2x + u2y) + 1e-100;
                
                curvatura = ((u2x*uyy) - (2*(ux*(uy*uxy))) + (u2y*uxx))/(((u2x + u2y)^(3/2))+(1e-100));
                if modelOption == 1
                    curvLtda = (curvatura*g(i,j)) + (1/mag)*((ux*gx(i,j)) + (uy*gy(i,j)));
                end
                
                delta(i,j) = 1;%(1/pi)*(EPSILON/(EPSILON^2 + u(i,j)^2));
                
                F = (delta(i,j).* ((VU*curvLtda) - LAMBDA*competition(i,j)) +1e-10);
                F = F / sqrt(sum(sum(F.^2)));%norm by square root of sum of square elements
                un(i,j) = u(i,j) - dt*F;
                
            end
        end
        y = double(u == 0);
        y2 = double(un == 0);
        residue = sum(sum(y.*y2));
        u = un;


      
         fprintf('============== ITERACAO %d ============\n\n\n',t);
%          if t == 1 || mod(t,10) == 0
%              %pathw = sprintf('%s_%s_M_%d.png',path,id,t);
%              m1 = mean(M1);
%              m2 = mean(M2);
%              temp = m1* HEps(u,EPSILON) + m2*(1- HEps(u,EPSILON));
%              %imwrite(uint8(temp),pathw);
%              pathw = sprintf('%slevelset_%d.png',outputPath,t);
%              level = (m1+m2)/2;
% %              figure,imshow(u,[]);
% %              m1
% %              m2
% %              level
% %              pause;
%              tempBW = temp > level;
%              bw = edge(tempBW ,'canny');
%              if canais == 3
%                  imgOut = Iout;
%              else
%                  imgOut = I;
%              end
%              imgBoundaries = zeros(altura,largura);
%              for i = 14:altura-13
%                 for j= 14:largura-13
%                     if bw(i,j) == 1;
%                         if canais == 3
%                             imgOut(i,j,1) = 255;
%                             imgOut(i,j,2) = 255;
%                             imgOut(i,j,3) = 255;
%                         else
%                             imgOut(i,j) = 255;
%                         end
%                         imgBoundaries(i,j) = 1;
%                     end
%                 end
%              end
%              imwrite(uint8(imgOut),pathw);
%              
%          end
    end
    
    m1 = mean(M1);
    m2 = mean(M2);
    temp = m1* HEps(u,EPSILON) + m2*(1- HEps(u,EPSILON));
    
    %pathw = sprintf('%slevelset_%d.png',outputPath,t);
    level = (m1+m2)/2;
    tempBW = temp > level;
    bw = edge(tempBW ,'canny');
     imgBoundaries = zeros(altura,largura);
             for i = 14:altura-13
                for j= 14:largura-13
                    if bw(i,j) == 1
                        if canais == 3
                            imgOut(i,j,1) = 255;
                            imgOut(i,j,2) = 255;
                            imgOut(i,j,3) = 255;
                        else
                            imgOut(i,j) = 255;
                        end
                        imgBoundaries(i,j) = 1;
                    end
                end
             end
    %pathw = sprintf('%ssegmentacao.png',);
    imwrite(imgBoundaries,outputPath);
   
    imgBoundaries2 = imfill(imgBoundaries, 'holes');
    
    todasBordas = bwboundaries(imgBoundaries2,'noholes');
    length(todasBordas)

    for k = 1:length(todasBordas)
        
        x = todasBordas{k}(:,1);
        %y = todasBordas{k}(:,2);
        
        tams = length(x);
        
%         fronteiras = zeros(altura,largura);
%         
%         for l = 1:tamx
%             fronteiras(x(l),y(l)) = 1;
%         end
    end
    
    [ordenado,indices] = sort(tams,'descend');
    
    %tams(k) = tamx;
    
    x = todasBordas{indices(1)}(:,1);
    y = todasBordas{indices(1)}(:,2);
    tamx = length(x);
    
    fronteiras = zeros(altura,largura);
    for l = 1:tamx
        fronteiras(x(l),y(l)) = 1;
    end
    
    bola = strel('disk',2);
    bin = imdilate(imgBoundaries,bola);
    fronteiras = imerode(bin,bola);
    
    imgFilled = imfill(fronteiras, 'holes');

    % Preenche a regiao internamente ao contorno
    area = imfill(imgFilled, 'holes');
    if tamx < 100
        area = imerode(area,bola);
    end
    %areaR = imresize(area,2);
    %area = areaR > 0.2;

    pathw = sprintf('%s%s.png',segPath,imgName);
    imwrite(area,pathw);
    %pathBin = sprintf('binary/%s_bin.png',lista(i).name(1:end-4));
    %imwrite(area,pathw);

end

function m = verificaLimites(x,lim)
    
    tamanho = length(x);
        
    m = 0;
    for i = 1:tamanho
        if x(i) == 1 | x(i) == lim
            m = i;
        end
    end
end


function mind = distCurva(x,y,curva)
    %(x,y) : coordenadas do ponto P
    %Vetor de pontos da curva
    
    for i = 1:length(curva)
        d(i) = sqrt((x - curva(i,1))^2 + (y - curva(i,2))^2);
    end
    
    mind = min(d);
end

function out = signLSM(x)
    if x > 0
        out = 1;
    else
        if x < 0
            out = -1;
        else
            out = 0;
        end
    end
end

function dprob = p(Img,m,DESVIO)

    pi = 3.1416;
    
    [altura,largura,canais]=size(Img);
    
    firstTerm = 1 / (sqrt(2*pi)*DESVIO);
    eTerm = 2*(DESVIO^2);
    
    for i = 1:altura
        for j = 1:largura

            x = Img(i,j);

            dprob(i,j) = firstTerm * exp( - (  (x-m)^2 ) / eTerm ) ;
            
        end
    end
end

function out = HEps(x,EPSILON)
    out = 0.5*(1+((2/pi)*atan(x/EPSILON)));
end

function prob = multiGaussian(x,M,Cov)
    
    %Define pi
    pi = 3.1416;

    %Calcula a dimensao do espaco
    d = length(M);
    
    detCov = det(Cov);
    
    baixo = ((2*pi)^(d/2))*(detCov^0.5);
    frac1 = 1/baixo;
    
    aT = x - M;
    a = (x - M)';
    invCov = inv(Cov);
    
    partexp = -0.5*aT*invCov*a;
    
    prob = frac1*partexp;

end

% function F = propagacao(u)
%     [a l] = size(u);
%     F = zeros(a,l);
%     menor = min(min(u));
%     maior = max(max(u));
%     
%     for j =2:a-1
%         for i = 2:l-1
%             if u(j,i) <= 60
%                 F(j,i) = u(j,i) - menor;
%             else
%                 F(j,i) = maior - u(j,i);
%             end
%         end
%     end
%           
% end
