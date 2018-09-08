function [z1,z2,epsilon,delta] = FactorScoresEst(Data,NumbObsVar)
%%function [z1,z2] = FactorScoresEst(Data,NumbObsVar)
  
  %{
With this function you can compute the Factor Scores of a latent variable
model. 

Input
DATA: This is one file with all observable variables for the independent
and the dependent variables. The observable variables for the dependent
variable have to be the last ones.
How it should be structured can be seen in the file "TestData" where 
datasets are simulated which you can use.
NUMBOBSVAR: This is a vector giving the number of observable variables for
every latent variable in the order they are given in the dataset. For
example [3,4,3] if there are 2 independent latent variables with 3 and 4
observable variables respectively and if the dependent latent variable has 
3 observable variables.

Output
Z1: Is a cell array containing the estimated factor scores of the 
independent latent variable
Z2: Is a vector containing the estimated factor scores of the dependent
variable
EPSILON: Is the estimation of the errors in the measurement model for the
independent variables.
DELTA: Is the estimation of the errors in the measurement model for the
dependent variable.
  %}

	 %AnzahlMessungen ist im Paper n und stimmt in X und Y überein
  dimY = NumbObsVar(end);

				% Load Data
  Daten = dlmread(Data);



  %% %Bring Data in the right form (rows = manifest Var; column = person)
  DimUnVar = length(NumbObsVar)-1;
  X = cell(1,DimUnVar);
  NumbObsVar_1 = [0,NumbObsVar];
  if size(Daten,2) == sum(NumbObsVar)
    for i=1:DimUnVar
      X(i) = {Daten(:,((i-1)*sum(NumbObsVar_1(1:i))+1):(sum(NumbObsVar_1(1:(i+1)))))};
    end
    Y = Daten(:,((end-NumbObsVar(end)+1):end));
    X = cellfun(@transpose,X,'UniformOutput',false);
    Y = Y';
  else
    for i=1:DimUnVar
      X(i) = {Daten((sum(NumbObsVar_1(1:i))+1):(sum(NumbObsVar_1(1:(i+1)))),:)};
    end
    Y = Daten(((end-NumbObsVar(end)+1):end),:);
  end
  clear i
  NumbObs = size(Y,2);

  %% %Conditions fulfilled of the Paper?
  for it = 1:DimUnVar
    if NumbObsVar(it) < 3
      disp('Dimension of one latent variable is too small (must be 3 or higher)');
      return
    end
  end
  if dimY < 3
    disp('Dimension of Y is too small (must be 3 or higher)');
    return
  end

  %% Estimation of coefficients of measurement model
  a = zeros(sum(NumbObsVar(1:(end-1))),1);
  for it = 1:DimUnVar
    a(sum(NumbObsVar_1(1:it))+1,1) = 1;
    a(sum(NumbObsVar_1(1:it))+2,1) = ((1/NumbObsVar_1(it+1))*sum(X{it}(2,:)...
								 *X{it}(3,:)'))/((1/NumbObsVar_1(it+1))*sum(X{it}(1,:)*X{it}(3,:)'));
    for indexa=3:NumbObsVar(it)
      a(sum(NumbObsVar_1(1:it))+indexa,1) = ((1/NumbObsVar_1(it+1))*...
					     sum(X{it}(2,:)*X{it}(indexa,:)'))/((1/NumbObsVar_1(it+1))*sum(X{it}(1,:)*X{it}(2,:)'));
    end
  end
  b = zeros(dimY,1);
  b(1,1) = 1;
  b(2,1) = ((1/dimY)*sum(Y(2,:)*Y(3,:)'))/((1/dimY)*sum(Y(1,:)*Y(3,:)'));
  for indexb=3:dimY
    b(indexb,1) = ((1/dimY)*sum(Y(2,:)*Y(indexb,:)'))/((1/dimY)*sum(Y(1,:)*Y(2,:)'));
  end
  %% % Estimation of Factor Scores
  
  Nn = round(NumbObs^(1/3));
  rho(1:Nn,1) = 1/Nn;
  alpha(1:Nn,1:(DimUnVar+1)) = randn(Nn,(DimUnVar+1));
  beta = randn(Nn,sum(NumbObsVar(1:(end-1))));
  gamma = randn(Nn,dimY);

  x0 = zeros((DimUnVar+1)*NumbObs,1);
  for it = 1:DimUnVar
    x0((it-1)*NumbObs+1:(it)*NumbObs,1) = (X{it}'*...
					   (1./a(sum(NumbObsVar_1(1:it))+1:sum(NumbObsVar_1(1:(it+1))),1)))/NumbObsVar(it);
  end
  x0(DimUnVar*NumbObs+1:(DimUnVar+1)*NumbObs) = (Y'*(1./b))/NumbObsVar(end);

  options = optimset('Algorithm','interior-point','MaxFunEvals',4500,'Display','off');

  [z,fval] = fmincon(@myfun,x0,[],[],[],[],[],[],@mycon,options);

  %% Definition of constraints
  function [c,ceq]=mycon(x)
    for it1 = 1:DimUnVar
      c(it) = (1/NumbObs)*sum(((x((it-1)*NumbObs+1:it*NumbObs)).^2))...
              - 1 - (1/NumbObs)*sum(((X{it}(1,:)).^2));
    end
    c(it+1) = (1/NumbObs)*sum(((x((it*NumbObs+1):...
				  (it+1)*NumbObs)).^2)) - 1 - ...
              (1/NumbObs)*sum(((Y(1,:)).^2));
    c(it+2) = (1/NumbObs)*sum(((x((it*NumbObs+1):...
				  (it+1)*NumbObs)).^4)) - 1 - (64/NumbObs)*...
							      sum(((Y(1,:)).^4)) - 72*(((1/NumbObs)*sum(Y(1,:)))^4);
    ceq = [];        
  end
  %% % Definition of objective function

  function Zielfunktionswert = myfun(x)
    Teil1 = 0;
    MatrixAb = zeros(NumbObs,sum(NumbObsVar_1));
    MatrixAc = zeros(NumbObs,sum(dimY));
    for r = 1:Nn
      TeilA = 0;
      TeilAa = ones(NumbObs,1);
      for i = 1:NumbObs
        for it2 = 1:DimUnVar+1
          TeilAa(i) = TeilAa(i)*1/(1+exp(NumbObs*(x(NumbObs*(it2-1)+i,1)-alpha(r,it2))));
        end
        for j1 = 1:DimUnVar
          for j2 = 1:NumbObsVar(j1)
            MatrixAb(i,sum(NumbObsVar_1(1:j1))+j2) = 1/(1+exp(NumbObs*((X{j1}(j2,i)-...
									a(sum(NumbObsVar_1(1:j1))+j2,1)*x(NumbObs*(j1-1)+i,1))-beta(r,sum(NumbObsVar_1(1:j1))+j2)))); 
          end
        end
        for k = 1:dimY % #ok<ALIGN>
          MatrixAc(i,k) = 1/(1+exp(NumbObs*((Y(k,i)-b(k,1)*x(end-NumbObs+i,1))-gamma(r,k))));
        end
        TeilAb = prod(MatrixAb(i,:));
        TeilAc = prod(MatrixAc(i,:));
        TeilA = TeilA + TeilAa(i) * TeilAb * TeilAc;
      end
      
      TeilA = TeilA/NumbObs;
      TeilBb = 1;
      TeilBc = 1;
      
      for j = 1:sum(NumbObsVar(1:(end-1))) % #ok<ALIGN>
	TeilBb = TeilBb * (1/NumbObs) * sum(MatrixAb(:,j)); 
      end
      for k = 1:dimY % #ok<ALIGN>
	TeilBc = TeilBc * (1/NumbObs) * sum(MatrixAc(:,k));
      end
      
      TeilBa = sum(TeilAa)/NumbObs;
      TeilB = TeilBa * TeilBb * TeilBc;
      Teil1a = (abs(TeilA-TeilB))^2;
      Teil1 = Teil1 + Teil1a*rho(r);

      %% clear MatrixAb MatrixAc
    end
    
    Teil2 = 0;
    Teil3 = 0;

    for j1 = 1:DimUnVar
      for j2 = 1:NumbObsVar(j1)
	Teil2a = 0;
	for i = 1:NumbObs
          Teil2a = Teil2a + (X{j1}(j2,i)-a(sum(NumbObsVar_1(1:j1))+j2,1)*x(NumbObs*(j1-1)+i,1));
	end
	Teil2 = Teil2 + (Teil2a*(1/NumbObs))^2;
      end
    end

    for k = 1:dimY
      Teil3a = 0;
      for i = 1:NumbObs
        Teil3a = Teil3a + (Y(k,i) - b(k,1)*x(end-NumbObs + i,1));
      end
      Teil3 = Teil3 + (Teil3a*(1/NumbObs))^2;
    end

    Teil = Teil1 + Teil2 + Teil3;


    Zielfunktionswert = Teil;
  end


  clear alpha beta gamma
  Zielfunktionswerte(1,1) = fval;

  %% % Make sure that the found Minimum is not local
% Here I compute the difference between the z-values. If they do no not
% change anymore, I assume that I found the global minimum

  z_all(:,1) = z;

      % with this criterion z may only differ after the second decimal
  Kriterium = 1/10000;
  Veraenderung = 10;
  NumbIt = 2;

  while NumbIt <= 10 && Veraenderung > Kriterium        
    
    alpha(1:Nn,1:(DimUnVar+1)) = randn(Nn,(DimUnVar+1));
    beta = randn(Nn,sum(NumbObsVar(1:(end-1))));
    gamma = randn(Nn,dimY);
    [z,fval] = fmincon(@myfun,z,[],[],[],[],[],[],@mycon,options);
    
    Zielfunktionswerte(NumbIt,1) = fval;
    z_all(:,NumbIt) = z;
    
    Veraenderung = sum((z_all(:,NumbIt)-z_all(:,NumbIt-1)).^2)/length(z);
    
    NumbIt = NumbIt + 1;
  end
  
  [~,min_Ind] = min(Zielfunktionswerte);

  z1 = cell(1,DimUnVar);
  for it3 = 1:DimUnVar
    z1{it3} = z_all((NumbObs*(it3-1)+1):(NumbObs*it3),min_Ind);
  end
  z2 = z_all((end-NumbObs+1):end,min_Ind);


  %% %estimation of the error
  for ind = 1:DimUnVar
    for indexj = 1:NumbObsVar(ind)
      for indexi = 1:NumbObs
        epsilon(sum(NumbObsVar_1(1:ind))+indexj, indexi) = X{ind}(indexj,indexi) - a(sum(NumbObsVar_1(1:ind))+indexj,1)*z1{ind}(indexi,1);
      end
    end
  end
  for indexk = 1:dimY
    for indexi = 1:NumbObs
      delta(indexk, indexi) = Y(indexk,indexi)-b(indexk,1)*z2(indexi,1);
    end
  end
  
end



