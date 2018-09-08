function [Factores,splinedeg,Knopt,Knots,z1] = BSplineEst(z1_inp,z2,DimUnVar)
 
%{
Input
Z1: Is a vector containing the estimated factor scores of the 
independent latent variable.
Z2: Is a vector containing the estimated factor scores of the dependent
variable
DimUnVar: This is a vector giving the number of observable variables for
every latent variable in the order they are given in the dataset. For
example [3,4,3] if there are 2 independent latent variables with 3 and 4
observable variables respectively and if the dependent latent variable has 
3 observable variables.

Output
FACTORES: Are the factores of the B-Splines
SPLINEDEG: Is the degree of the B-Splines
KNOTS: Are the used knots of the B-Splines
KNOPT: Is the optimal Number of Knots

With the results you can still use the PLOTS of the function ModelEst

%}

NumbObs = length(z2);

if iscell(z1_inp) == 0
   if size(z1_inp,2) == DimUnVar
       z1 = cell(1,DimUnVar);
       for it3 = 1:DimUnVar
           z1{it3} = z1_inp(:,it3);
       end
   elseif size(z1_inp,1) == DimUnVar
       z1 = cell(1,DimUnVar);
       for it3 = 1:DimUnVar
           z1{it3} = z1_inp(it3,:)';
       end
   elseif size(z1_inp,1) == length(z2)*DimUnVar
       z1 = cell(1,DimUnVar);
       for it3 = 1:DimUnVar
           z1{it3} = z1_inp((NumbObs*(it3-1)+1):(NumbObs*it3),1);
       end
   elseif size(z1_inp,2) == length(z2)*DimUnVar
       z1 = cell(1,DimUnVar);
       for it3 = 1:DimUnVar
           z1{it3} = z1_inp(1,(NumbObs*(it3-1)+1):(NumbObs*it3))';
       end
   end
else
    z1 = z1_inp;
end

% z1 = cell(1,DimUnVar);
% for it3 = 1:DimUnVar
%     z1{it3} = z_all((NumbObs*(it3-1)+1):(NumbObs*it3),min_Ind);
% end
% z2 = z_all((end-NumbObs+1):end,min_Ind);



if DimUnVar == 1
    z1 = z1{1};
    zmin = min(z1);
    zmax = max(z1);

    Trainingz1 = z1(1:round(length(z1)/2));
    Evaluationz1 = z1(round(length(z1)/2+1:end));
    Trainingz2 = z2(1:round(length(z2)/2));
    Evaluationz2 = z2(round(length(z2)/2)+1:length(z2));
    Grad = 4;

    Kn(1) = 0;

    for n = 1:11
        Kn(n+1) = round(n+Grad-2);
        if Kn(n+1) ~= Kn(n)
            diff = (zmax - zmin)/Kn(n+1);
            t = (zmin-(Grad+1)*diff):diff:(zmax+(Grad+1)*diff);
            B = bspline_basismatrix(Grad,t,Trainingz1);

            C = B'*B;
            d = B'*Trainingz2;
            C = sparse(C);
            [a,~] = lsqr(C,d,[],20);

            clear B b
        
            B = bspline_basismatrix(Grad,t,Evaluationz1);
            Geschaetzt = B*a;
            Geschaetzt = sum(Geschaetzt,2);
    
            Fehler = (Evaluationz2 - Geschaetzt).^2;
            Fehler = sum(Fehler)/length(z1);
    
            if n == 1
                kleinsterFehler = Fehler; 
                Factores = a;
                Knots = t;
                splinedeg = Grad-1;
                Knopt = Kn(n+1);
            end
    
            if Fehler < kleinsterFehler
                kleinsterFehler = Fehler;
                Factores = a;
                Knots = t;
                splinedeg = Grad-1;
                Knopt = Kn(n+1);
            end
    
            clear Fehler a C d B Fehler2 a2 C2 d2 b
        end
    end

    z1_t = cell(1,1);
    z1_t{1} = z1;
    clear z1
    z1 = z1_t;

elseif DimUnVar >= 2
    zmin = zeros(DimUnVar,1);
    zmax = zeros(DimUnVar,1);
    for i = 1:DimUnVar
       zmin(i) = min(z1{i}); 
       zmax(i) = max(z1{i}); 
    end

    Trainingz1 = cell(1,DimUnVar);
    Evaluationz1 = cell(1,DimUnVar);
    for it = 1:DimUnVar
        Trainingz1{it} = z1{it}(1:round(length(z1{1})/2));
        Evaluationz1{it} = z1{it}(round(length(z1{1})/2+1:end));
    end
    Trainingz2 = z2(1:round(length(z2)/2));
    Evaluationz2 = z2(round(length(z2)/2)+1:length(z2));

    Grad = 4;
    Kn = zeros(11,1);
    diff = zeros(DimUnVar,1);
    for n = 1:11
        Kn(n) = n+Grad-2;
        t = cell(1,DimUnVar);
        for i = 1:DimUnVar
            diff(i) = (zmax(i) - zmin(i))/Kn(n);        
            t{i} = (zmin(i)-(Grad+1)*diff(i)):diff(i):(zmax(i)+(Grad+1)*diff(i));
        end            
        B = [];
        B_ind = cell(1,DimUnVar);
        for i = 1:DimUnVar
            B_ind{i} = bspline_basismatrix(Grad,t{i},Trainingz1{i});
        end
        
        for it = 1:length(Trainingz1{1})
            B_it = B_ind{1}(it,:);
            for i = 2:DimUnVar
                B_it = B_it'*B_ind{i}(it,:);
                B_it = reshape(B_it,1,[]);
            end
            B = [B;B_it];
        end
        
        C = B'*B;
        d = B'*Trainingz2;
        C = sparse(C);
        [a,~] = lsqr(C,d,[],20);

        clear B b

        B = [];
        B_ind = cell(1,DimUnVar);
        for i = 1:DimUnVar
            B_ind{i} = bspline_basismatrix(Grad,t{i},Evaluationz1{i});
        end
        for it = 1:length(Evaluationz1{1})
            B_it = B_ind{1}(it,:);
            for i = 2:DimUnVar
                B_it = B_it'*B_ind{i}(it,:);
                B_it = reshape(B_it,1,[]);
            end
            B = [B;B_it];
        end

        Geschaetzt = B*a;
        Geschaetzt = sum(Geschaetzt,2);

        Fehler = (Evaluationz2 - Geschaetzt).^2;
        Fehler = sum(Fehler)/length(z2);

        if n == 1
           kleinsterFehler = Fehler; 
           Factores = a;
           Knots = t;
           splinedeg = Grad-1;
           Knopt = Kn(n);
        end

        if Fehler < kleinsterFehler
           kleinsterFehler = Fehler;
           Factores = a;
           Knots = t;
           splinedeg = Grad-1;
           Knopt = Kn(n);
        end

        clear Fehler a C d B Fehler2 a2 C2 d2 b
        
    end
end
