% Algorithme codé à partir de l'article 'PatchMatch: A Randomized Correspondence Algorithm for Structural Image Editing'
%   de Connelly Barnes, Eli Shechtman, Adam Finkelstein, Dan B Goldman
%   Princeton University, Adobe Systems, University of Washington

function [pp_voisins, distances] = PatchMatch_article(img, taille_patch)
% ARGUMENTS:
%   img: image en entrée
%   taille_patch: taille des patchs à rechercher, on choisit une taille de
%       9 par défaut
% OUTPUT:
%   pp_voisins: plus proches voisins sous la forme d'indices
%   distances: distance entre les patchs correpondant en norme L^2

%% Préparation des images en amont des opérations
% Taille de l'image en entrée
taille_img = size(img);

% On vérifie que la taille du patch est bien impaire sinon, on sort de la
%   fonction
if mod(taille_patch,2)==1
    largeur = (taille_patch - 1)/2;
else
    disp('Attention: la taille du patch doit être impaire pour faire fonctionner lalgorithme')
    return
end

% Conversion au format 'double' pour ne pas avoir de problème dans la
%   manipulation de calculs
img = double(img);

%% Initialisation
% Paramètres de la boucle itérative
max_iterations = 5; % Nombre maximal d'itérations

% Initialisation aléatoire des correspondances entre patchs
pp_voisins = cat(3, randi([1,taille_img(1)],taille_img(1:2)),...
    randi([1,taille_img(2)],taille_img(1:2)));
% Calculs des distances entre patchs associés
distances = inf(taille_img(1),taille_img(2));
for i = 1:taille_img(1)
    for j = 1:taille_img(2)
        int = -largeur:largeur;
        test_1 = (i+int>=1) & (i+int<=taille_img(1)) & (pp_voisins(i,j,1)+int>=1) & (pp_voisins(i,j,1)+int<=taille_img(1));
        test_2 = (j+int>=1) & (j+int<=taille_img(2)) & (pp_voisins(i,j,2)+int>=1) & (pp_voisins(i,j,2)+int<=taille_img(2));
        
        distance_patchs_rgb = img(i+int(test_1),j+int(test_2),:)...
            - img(pp_voisins(i,j,1)+int(test_1),pp_voisins(i,j,2)+int(test_2),:);
        distance_patchs_rgb = distance_patchs_rgb(:);
        penalite_rgb = sum(distance_patchs_rgb.^2)/length(3*distance_patchs_rgb);
        
        distances(i,j) = penalite_rgb;
    end
end

%% Boucle de calculs : Propagation puis Recherche Aléatoire
for iteration = 1:max_iterations
    disp(['Itération numéro: ',num2str(iteration)])
    
    % L'étape de propagation dépend de la parité de l'itération
    parite = mod(iteration,2)==1;
    if parite == 1
        i_seq = 1:taille_img(1);
        j_seq = 1:taille_img(2);
    else
        i_seq = taille_img(1):-1:1;
        j_seq = taille_img(2):-1:1;
    end
    
    for i = i_seq
      for j = j_seq
        %% Propagation
        if parite
            %%% Centre, Bas, Gauche
            prop_distances = inf(5,1);
            % Centre
            prop_distances(1) = distances(i,j);
            % Gauche
            if i-1>0 && pp_voisins(i-1,j,1)-1>0
                int = -largeur:largeur;
                test_1 = (i+int>=1) & (i+int<=taille_img(1)) & (pp_voisins(i-1,j,1)-1+int>=1) & (pp_voisins(i-1,j,1)-1+int<=taille_img(1));
                test_2 = (j+int>=1) & (j+int<=taille_img(2)) & (pp_voisins(i-1,j,2)+int>=1) & (pp_voisins(i-1,j,2)+int<=taille_img(2));

                distance_patchs_rgb = img(i+int(test_1),j+int(test_2),:)...
                    - img(pp_voisins(i-1,j,1)-1+int(test_1),pp_voisins(i-1,j,2)+int(test_2),:);
                distance_patchs_rgb = distance_patchs_rgb(:);
                penalite_rgb = sum(distance_patchs_rgb.^2)/length(3*distance_patchs_rgb);
        
                prop_distances(2) = penalite_rgb;
            end
            % Haut
            if j-1>0 && pp_voisins(i,j-1,2)-1>0
                int = -largeur:largeur;
                test_1 = (i+int>=1) & (i+int<=taille_img(1)) & (pp_voisins(i,j-1,1)+int>=1) & (pp_voisins(i,j-1,1)+int<=taille_img(1));
                test_2 = (j+int>=1) & (j+int<=taille_img(2)) & (pp_voisins(i,j-1,2)-1+int>=1) & (pp_voisins(i,j-1,2)-1+int<=taille_img(2));

                distance_patchs_rgb = img(i+int(test_1),j+int(test_2),:)...
                    - img(pp_voisins(i,j-1,1)+int(test_1),pp_voisins(i,j-1,2)-1+int(test_2),:);
                distance_patchs_rgb = distance_patchs_rgb(:);
                penalite_rgb = sum(distance_patchs_rgb.^2)/length(3*distance_patchs_rgb);
        
                prop_distances(3) = penalite_rgb;
            end
            % Rotation Gauche
            if i-2>0
                f = 2*[pp_voisins(i-1,j,1)-(i-1),pp_voisins(i-1,j,2)-j]-[pp_voisins(i-2,j,1)-(i-2),pp_voisins(i-2,j,2)-j];
                if i+f(1)>0 && i+f(1)<=taille_img(1) && j+f(2)>0 && j+f(2)<=taille_img(2)
                    int = -largeur:largeur;
                    test_1 = (i+int>=1) & (i+int<=taille_img(1)) & (i+f(1)+int>=1) & (i+f(1)+int<=taille_img(1));
                    test_2 = (j+int>=1) & (j+int<=taille_img(2)) & (j+f(2)+int>=1) & (j+f(2)+int<=taille_img(2));

                    distance_patchs_rgb = img(i+int(test_1),j+int(test_2),:)...
                        - img(i+f(1)+int(test_1),j+f(2)+int(test_2),:);
                    distance_patchs_rgb = distance_patchs_rgb(:);
                    penalite_rgb = sum(distance_patchs_rgb.^2)/length(3*distance_patchs_rgb);

                    prop_distances(4) = penalite_rgb;
                end
            end
            % Rotation Haut
            if j-2>0
                f = 2*[pp_voisins(i,j-1,1)-i,pp_voisins(i,j-1,2)-(j-1)]-[pp_voisins(i,j-2,1)-i,pp_voisins(i,j-2,2)-(j-2)];
                if i+f(1)>0 && i+f(1)<=taille_img(1) && j+f(2)>0 && j+f(2)<=taille_img(2)
                    int = -largeur:largeur;
                    test_1 = (i+int>=1) & (i+int<=taille_img(1)) & (i+f(1)+int>=1) & (i+f(1)+int<=taille_img(1));
                    test_2 = (j+int>=1) & (j+int<=taille_img(2)) & (j+f(2)+int>=1) & (j+f(2)+int<=taille_img(2));

                    distance_patchs_rgb = img(i+int(test_1),j+int(test_2),:)...
                        - img(i+f(1)+int(test_1),j+f(2)+int(test_2),:);
                    distance_patchs_rgb = distance_patchs_rgb(:);
                    penalite_rgb = sum(distance_patchs_rgb.^2)/length(3*distance_patchs_rgb);

                    prop_distances(5) = penalite_rgb;
                end
            end
            
            [~,ind] = min(prop_distances);
            % Propagation à Gauche
            if ind == 2
                % Remplacement du plus proche voisin
                pp_voisins(i,j,:) = pp_voisins(i-1,j,:);
                if pp_voisins(i,j,1)-1>0 && pp_voisins(i,j,1)-1<=taille_img(1)
                    pp_voisins(i,j,1) = pp_voisins(i,j,1)-1;
                end
                distances(i,j) = prop_distances(2);
            % Propagation en Haut
            elseif ind == 3
                % Remplacement du plus proche voisin
                pp_voisins(i,j,:) = pp_voisins(i,j-1,:);
                if pp_voisins(i,j,2)-1>0 && pp_voisins(i,j,2)-1<=taille_img(2)
                    pp_voisins(i,j,2) = pp_voisins(i,j,2)-1;
                end
                distances(i,j) = prop_distances(3);
            % Propagarion Rotation Gauche
            elseif ind == 4
                f = 2*[pp_voisins(i-1,j,1)-(i-1),pp_voisins(i-1,j,2)-j]-[pp_voisins(i-2,j,1)-(i-2),pp_voisins(i-2,j,2)-j];
                % Remplacement du plus proche voisin
                pp_voisins(i,j,1) = i+f(1);
                pp_voisins(i,j,2) = j+f(2);
                distances(i,j) = prop_distances(4);
            % Propagarion Rotation Haut
            elseif ind == 5
                f = 2*[pp_voisins(i,j-1,1)-i,pp_voisins(i,j-1,2)-(j-1)]-[pp_voisins(i,j-2,1)-i,pp_voisins(i,j-2,2)-(j-2)];
                % Remplacement du plus proche voisin
                pp_voisins(i,j,1) = i+f(1);
                pp_voisins(i,j,2) = j+f(2);
                distances(i,j) = prop_distances(5);
            end
            
        else
            %%% Centre, Droite, Bas
            prop_distances = inf(5,1);
            % Centre
            prop_distances(1) = distances(i,j);
            % Droite
            if i+1<taille_img(1) && pp_voisins(i+1,j,1)+1<taille_img(1)
                int = -largeur:largeur;
                test_1 = (i+int>=1) & (i+int<=taille_img(1)) & (pp_voisins(i+1,j,1)+1+int>=1) & (pp_voisins(i+1,j,1)+1+int<=taille_img(1));
                test_2 = (j+int>=1) & (j+int<=taille_img(2)) & (pp_voisins(i+1,j,2)+int>=1) & (pp_voisins(i+1,j,2)+int<=taille_img(2));

                distance_patchs_rgb = img(i+int(test_1),j+int(test_2),:)...
                    - img(pp_voisins(i+1,j,1)+1+int(test_1),pp_voisins(i+1,j,2)+int(test_2),:);
                distance_patchs_rgb = distance_patchs_rgb(:);
                penalite_rgb = sum(distance_patchs_rgb.^2)/length(3*distance_patchs_rgb);
        
                prop_distances(2) = penalite_rgb;
            end
            % Bas
            if j+1<taille_img(2) && pp_voisins(i,j+1,1)+1<taille_img(2)
                int = -largeur:largeur;
                test_1 = (i+int>=1) & (i+int<=taille_img(1)) & (pp_voisins(i,j+1,1)+int>=1) & (pp_voisins(i,j+1,1)+int<=taille_img(1));
                test_2 = (j+int>=1) & (j+int<=taille_img(2)) & (pp_voisins(i,j+1,2)+1+int>=1) & (pp_voisins(i,j+1,2)+1+int<=taille_img(2));

                distance_patchs_rgb = img(i+int(test_1),j+int(test_2),:)...
                    - img(pp_voisins(i,j+1,1)+int(test_1),pp_voisins(i,j+1,2)+1+int(test_2),:);
                distance_patchs_rgb = distance_patchs_rgb(:);
                penalite_rgb = sum(distance_patchs_rgb.^2)/length(3*distance_patchs_rgb);
        
                prop_distances(3) = penalite_rgb;
            end
            % Rotation Droite
            if i+2<=taille_img(1)
                f = 2*[pp_voisins(i+1,j,1)-(i+1),pp_voisins(i+1,j,2)-j]-[pp_voisins(i+2,j,1)-(i+2),pp_voisins(i+2,j,2)-j];
                if i+f(1)>0 && i+f(1)<=taille_img(1) && j+f(2)>0 && j+f(2)<=taille_img(2)
                    int = -largeur:largeur;
                    test_1 = (i+int>=1) & (i+int<=taille_img(1)) & (i+f(1)+int>=1) & (i+f(1)+int<=taille_img(1));
                    test_2 = (j+int>=1) & (j+int<=taille_img(2)) & (j+f(2)+int>=1) & (j+f(2)+int<=taille_img(2));

                    distance_patchs_rgb = img(i+int(test_1),j+int(test_2),:)...
                        - img(i+f(1)+int(test_1),j+f(2)+int(test_2),:);
                    distance_patchs_rgb = distance_patchs_rgb(:);
                    penalite_rgb = sum(distance_patchs_rgb.^2)/length(3*distance_patchs_rgb);

                    prop_distances(4) = penalite_rgb;
                end
            end
            % Rotation Bas
            if j+2<=taille_img(2)
                f = 2*[pp_voisins(i,j+1,1)-i,pp_voisins(i,j+1,2)-(j+1)]-[pp_voisins(i,j+2,1)-i,pp_voisins(i,j+2,2)-(j+2)];
                if i+f(1)>0 && i+f(1)<=taille_img(1) && j+f(2)>0 && j+f(2)<=taille_img(2)
                    int = -largeur:largeur;
                    test_1 = (i+int>=1) & (i+int<=taille_img(1)) & (i+f(1)+int>=1) & (i+f(1)+int<=taille_img(1));
                    test_2 = (j+int>=1) & (j+int<=taille_img(2)) & (j+f(2)+int>=1) & (j+f(2)+int<=taille_img(2));

                    distance_patchs_rgb = img(i+int(test_1),j+int(test_2),:)...
                        - img(i+f(1)+int(test_1),j+f(2)+int(test_2),:);
                    distance_patchs_rgb = distance_patchs_rgb(:);
                    penalite_rgb = sum(distance_patchs_rgb.^2)/length(3*distance_patchs_rgb);

                    prop_distances(5) = penalite_rgb;
                end
            end
            
            [~,ind] = min(prop_distances);
            % Propagation à Droite
            if ind == 2
                % Remplacement du plus proche voisin
                pp_voisins(i,j,:) = pp_voisins(i+1,j,:);
                if pp_voisins(i,j,1)+1>0 && pp_voisins(i,j,1)+1<=taille_img(1)
                    pp_voisins(i,j,1) = pp_voisins(i,j,1)+1;
                end
                distances(i,j) = prop_distances(2);
            % Propagation en Bas
            elseif ind == 3
                % Remplacement du plus proche voisin
                pp_voisins(i,j,:) = pp_voisins(i,j+1,:);
                if pp_voisins(i,j,2)+1>0 && pp_voisins(i,j,2)+1<=taille_img(2)
                    pp_voisins(i,j,2) = pp_voisins(i,j,2)+1;
                end
                distances(i,j) = prop_distances(3);
            % Propagation Rotation Droite
            elseif ind == 4
                f = 2*[pp_voisins(i+1,j,1)-(i+1),pp_voisins(i+1,j,2)-j]-[pp_voisins(i+2,j,1)-(i+2),pp_voisins(i+2,j,2)-j];
                % Remplacement du plus proche voisin
                pp_voisins(i,j,1) = i+f(1);
                pp_voisins(i,j,2) = j+f(2);
                distances(i,j) = prop_distances(4);
            % Propagation Rotation Bas
            elseif ind == 5
                f = 2*[pp_voisins(i,j+1,1)-i,pp_voisins(i,j+1,2)-(j+1)]-[pp_voisins(i,j+2,1)-i,pp_voisins(i,j+2,2)-(j+2)];
                % Remplacement du plus proche voisin
                pp_voisins(i,j,1) = i+f(1);
                pp_voisins(i,j,2) = j+f(2);
                distances(i,j) = prop_distances(5);
            end
        end

        %% Recherche aléatoire
        
        % Paramètres de la recherche aléatoire
        L = 10;
        w = max(size(img));
        alpha = 1/2;
        
        % Initialisation de la Recherche aléatoire
        R_i_x = randi([-1,1],L,1);
        R_i_y = randi([-1,1],L,1);
        distance = inf(1+L,1);
        distance(1) = distances(i,j);
        
        for it = 1:L
            if pp_voisins(i,j,1)+round(alpha^(it-1)*w*R_i_x(it))<=0 || pp_voisins(i,j,1)+round(alpha^(it-1)*w*R_i_x(it))>=taille_img(1) ...
                    || pp_voisins(i,j,2)+round(alpha^(it-1)*w*R_i_y(it))<=0 || pp_voisins(i,j,2)+round(alpha^(it-1)*w*R_i_y(it))>=taille_img(2)
                distance(it+1) = nan;
            else
                int = -largeur:largeur;
                test_1 = (i+int>=1) & (i+int<=taille_img(1)) & (pp_voisins(i,j,1)+int+round(alpha^(it-1)*w*R_i_x(it))>=1) & (pp_voisins(i,j,1)+int+round(alpha^(it-1)*w*R_i_x(it))<=taille_img(1));
                test_2 = (j+int>=1) & (j+int<=taille_img(2)) & (pp_voisins(i,j,2)+int+round(alpha^(it-1)*w*R_i_y(it))>=1) & (pp_voisins(i,j,2)+int+round(alpha^(it-1)*w*R_i_y(it))<=taille_img(2));

                distance_patchs_rgb = img(i+int(test_1),j+int(test_2),:)...
                    - img(pp_voisins(i,j,1)+int(test_1)+round(alpha^(it-1)*w*R_i_x(it)),pp_voisins(i,j,2)+int(test_2)+round(alpha^(it-1)*w*R_i_y(it)),:);
                distance_patchs_rgb = distance_patchs_rgb(:);
                penalite_rgb = sum(distance_patchs_rgb.^2)/length(3*distance_patchs_rgb);

                distance(it+1) = penalite_rgb;
            end
        end
        [~,ind] = sort(distance,'ascend');
        if length(distance(~isnan(distance)))<2
            ind = ind(1);
        else
            ind = ind(2);
        end
        distances(i,j) = distance(ind);
        if ind ~= 1
            pp_voisins(i,j,1) = pp_voisins(i,j,1) + round(alpha^(it-2)*w*R_i_x(it-1));
            pp_voisins(i,j,2) = pp_voisins(i,j,2)+round(alpha^(it-2)*w*R_i_y(it-1));
        end
      end
    end

end

end