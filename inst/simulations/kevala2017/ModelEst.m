function [z1,z2,epsilon,delta,factores,splinedeg,Knots] = ModelEst(Data,NumbObsVar,Plot)
  
%{
This is the main file you can use to estimate the factor scores and the
regression function of your given data. It uses the functions 
FactorScoresEst and BSplineEst.

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
PLOT: If you want to see the estimates plotted there are 3 options:
    1: You get a plot of the estimated factor scores as a function of each 
    other if there are not more than 3 latent variables
    2: You get a plot of the estimated factor scores of the independent
    latent variables as a function of the the dependent latent variable 
    estimated by the B-Splines (the regression function)
    3: You get both plots of options 1 and 2

Output
Z1: Is a cell array containing the estimated factor scores of the 
independent latent variable
Z2: Is a vector containing the estimated factor scores of the dependent
variable
FACTORES: Are the factores of the B-Splines
SPLINEDEG: Is the degree of the B-Splines
KNOTS: Are the used knots of the B-Splines
EPSILON: Is the estimation of the errors in the measurement model for the
independent variables.
DELTA: Is the estimation of the errors in the measurement model for the
dependent variable.
%}

%%

% Estimation of the Factor Scores
[z1,z2,epsilon,delta] = FactorScoresEst(Data,NumbObsVar);

% Estimation of the Regression between the Factor Scores
DimUnVar = length(NumbObsVar)-1;
[factores,splinedeg,Knopt,Knots,z1] = BSplineEst(z1,z2,DimUnVar);

%% % Possible Plots
% If you estimated only the regression with given Factor Scores you can
% execute just the lines of every case to get the plots
if nargin == 3
    switch Plot
        case 1  %Factor scores against each other
            if DimUnVar == 1
                figure('Name','IV_vs_DV','NumberTitle','off');
                plot(z1{1},z2,'o')
            elseif DimUnVar == 2
                figure('Name','IV_vs_DV','NumberTitle','off');
                plot3(z1{1},z1{2},z2,'o')
            else
                disp(['Your input data has more than 3 dimensions and can' ...
                    'not be plotted']);
                return
            end
        case 2  %independant variable against B-Spline estimates of dependant variable
            if DimUnVar == 1
                B = bspline_basismatrix(splinedeg+1,Knots,z1{1});
                Schaetzungen = B*factores;
                figure('Name','IV_vs_Estimates','NumberTitle','off');
                plot(z1{1},Schaetzungen,'o')
            elseif DimUnVar == 2
                B = [];
                B_ind = cell(1,DimUnVar);
                for it_B = 1:DimUnVar
                    B_ind{it_B} = bspline_basismatrix(splinedeg+1,Knots{it_B},z1{it_B});
                end
                for it = 1:length(z1{1})
                B_it = B_ind{1}(it,:);
                for it_B = 2:DimUnVar
                    B_it = B_it'*B_ind{it_B}(it,:);
                    B_it = reshape(B_it,1,[]);
                end
                B = [B;B_it];
                end
                Schaetzungen = B*factores;
                figure('Name','IV_vs_Estimates','NumberTitle','off');
                plot3(z1{1},z1{2},Schaetzungen,'o')
            else
                disp(['Your input data has more than 3 dimensions and can' ...
                    'not be plotted']);
                return
            end            
        case 3  %both plots
            if DimUnVar == 1
                B = bspline_basismatrix(splinedeg+1,Knots,z1{1});
                Schaetzungen = B*factores;
                figure('Name','IV_vs_Estimates','NumberTitle','off');
                plot(z1{1},Schaetzungen,'o')
                figure('Name','IV_vs_DV','NumberTitle','off');
                plot(z1{1},z2,'o')
            elseif DimUnVar == 2
                B = [];
                B_ind = cell(1,DimUnVar);
                for it_B = 1:DimUnVar
                    B_ind{it_B} = bspline_basismatrix(splinedeg+1,Knots{it_B},z1{it_B});
                end
                for it = 1:length(z1{1})
                B_it = B_ind{1}(it,:);
                for it_B = 2:DimUnVar
                    B_it = B_it'*B_ind{it_B}(it,:);
                    B_it = reshape(B_it,1,[]);
                end
                B = [B;B_it];
                end
                Schaetzungen = B*factores;
                figure('Name','IV_vs_Estimates','NumberTitle','off');
                plot3(z1{1},z1{2},Schaetzungen,'o')
                figure('Name','IV_vs_DV','NumberTitle','off');
                plot3(z1{1},z1{2},z2,'o')
            else
                disp(['Your input data has more than 3 dimensions and can' ...
                    'not be plotted']);
                return
            end            
    end
end
