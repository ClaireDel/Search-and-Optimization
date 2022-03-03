%% Assignment : Livestock feeding

clear all;
clc


% %% PART 1 and 2
% 
% %% Problem
% 
% % Cost function to minimize
% f = [0.30 0.90]; 
% 
% % Constraints matrix
% A = [-1     -1    ; 
%      -1.3   -1.4  ;
%       1.3    1.4  ;
%       0.21  -0.3  ;
%      -0.03   0.01 ;] ;
% 
% b = [-400 -600 1000 0 0]; 
%  
% %% Solving by Graphical Method
% 
% figure(1)
% lb = [0 0]; % lower bounds
% ub = [1000 1000]; % upper bound (for the plot)
% [val, fval] = linprog(f, A, b, [], [], lb, []);
% plotregion(-A,-b,lb,[],'g',0.5) % other file
% xlabel('x1: Corn'), ylabel('x2: Fodder')
% axis([lb(1) ub(1) lb(2) ub(2)]), grid
% title('Linear Problem for livestock feeding: Graph Solving')
% hold on
% 
% % Constraints lines:
% x1 = lb(1):1000; 
% x2 = lb(2):1000;
% x2_1 = 400 - x1; 
% x2_2 = 600/1.4 - 1.3/1.4*x1 ; 
% x2_3 = 1000/1.4 - 1.3/1.4*x1;
% x2_4 = 0.21/0.3*x1;
% x2_5 = 3*x1;
% 
% % Optimal cost function line:
% x2_cf = -(0.3.*x1) / 0.9; 
% obj = plot(x1, zeros(size(x1)), zeros(size(x2)), x2, x1, x2_1, x1, x2_2, x1, x2_3, x1, x2_4, x1, x2_5, x1, x2_cf);
% obj(3).LineStyle = '--';
% obj(4).LineStyle = '--';
% obj(5).LineStyle = '--';
% obj(6).LineStyle = '--';
% obj(7).LineStyle = '--';
% obj(8).LineStyle = '--';
% % plot optimal solution:
% plot(val(1), val(2), 'r*');
% hold off % releases graph
% % set axes and figure legend:
% legend(obj, {'x \geq 0', 'y \geq 0', ...
%     'c1: y \geq 400 - x', ...
%     'c2: y \geq 600/1.4 - 1.3/1.4x', ...
%     'c3: y \leq 1000/1.4 - 1.3/1.4x', ...
%     'c4: y \geq 0.21/0.3x', ...
%     'c5: y \leq 3x', ...
%     'cf: y = -0.3x / 0.9'}, ...
%     'Location', 'Best');
% 
% 
% %% Solving by linprog
% 
% [val, fval] = linprog(f, A, b, [], [], lb, []);
% disp('val = ');
% disp(val);
% disp('fval = ');
% disp(fval);
% % 
% 
% %% Sensitivity Analysis 1
% 
% % Cost function to minimize
% f = [0.30 0.90]; 
% 
% % New constraints matrix
% Abis = [-1     -1    ; 
%         -1.3   -1.4  ;
%          1.3    1.4  ;
%          0.15  -0.3  ;
%         -0.03   0.01 ;] ;
% 
% bbis = [-400 -600 1000 0 0]; 
% 
% figure(2)
% lb = [0 0]; % lower bounds
% ub = [1000 1000]; % upper bound (for the plot)
% [valbis, fvalbis] = linprog(f, Abis, bbis, [], [], lb, []);
% plotregion(-Abis,-bbis,lb,[],'g',0.5) % other file
% xlabel('x1: Corn'), ylabel('x2: Fodder')
% axis([lb(1) ub(1) lb(2) ub(2)]), grid
% title('Linear Problem for livestock feeding: Graph Solving')
% hold on
% 
% % Constraints lines:
% x1 = lb(1):1000; 
% x2 = lb(2):1000;
% x2_1 = 400 - x1; 
% x2_2 = 600/1.4 - 1.3/1.4*x1 ; 
% x2_3 = 1000/1.4 - 1.3/1.4*x1;
% x2_4 = 0.15/0.3*x1;
% x2_5 = 3*x1;
% 
% % Optimal cost function line:
% x2_cf = -(0.3.*x1) / 0.9; 
% obj = plot(x1, zeros(size(x1)), zeros(size(x2)), x2, x1, x2_1, x1, x2_2, x1, x2_3, x1, x2_4, x1, x2_5, x1, x2_cf);
% obj(3).LineStyle = '--';
% obj(4).LineStyle = '--';
% obj(5).LineStyle = '--';
% obj(6).LineStyle = '--';
% obj(7).LineStyle = '--';
% obj(8).LineStyle = '--';
% % plot optimal solution:
% plot(valbis(1), valbis(2), 'r*');
% hold off % releases graph
% % set axes and figure legend:
% legend(obj, {'x \geq 0', 'y \geq 0', ...
%     'c1: y \geq 400 - x', ...
%     'c2: y \geq 600/1.4 - 1.3/1.4 x', ...
%     'c3: y \leq 1000/1.4 - 1.3/1.4 x', ...
%     'c4: y \geq 0.15/0.3x', ...
%     'c5: y \leq 3x', ...
%     'cf: y = -0.3x / 0.9'}, ...
%     'Location', 'Best');
% 
% [valbis, fvalbis] = linprog(f, Abis, bbis, [], [], lb, []);
% disp('val = ');
% disp(valbis);
% disp('fval = ');
% disp(fvalbis);
% 
% %% Sensitivity Analysis 2
% 
% % Cost function to minimize
% f = [0.30 0.90]; 
% 
% % New constraints matrix
% Abis2 = [-1     -1    ; 
%         -1.3   -2  ;
%          1.3    2  ;
%          0.15  -0.3  ;
%         -0.03   0.01 ;] ;
% 
% bbis2 = [-400 -600 1000 0 0]; 
% 
% figure(3)
% lb = [0 0]; % lower bounds
% ub = [1000 1000]; % upper bound (for the plot)
% [valbis2, fvalbis2] = linprog(f, Abis2, bbis2, [], [], lb, []);
% plotregion(-Abis2,-bbis2,lb,[],'g',0.5) % other file
% xlabel('x1: Corn'), ylabel('x2: Fodder')
% axis([lb(1) ub(1) lb(2) ub(2)]), grid
% title('Linear Problem for livestock feeding: Graph Solving')
% hold on
% 
% % Constraints lines:
% x1 = lb(1):1000; 
% x2 = lb(2):1000;
% x2_1 = 400 - x1; 
% x2_2 = 600/2 - 1.3/2*x1 ; 
% x2_3 = 1000/2 - 1.3/2*x1;
% x2_4 = 0.15/0.3*x1;
% x2_5 = 3*x1;
% 
% % Optimal cost function line:
% x2_cf = -(0.3.*x1) / 0.9; 
% obj = plot(x1, zeros(size(x1)), zeros(size(x2)), x2, x1, x2_1, x1, x2_2, x1, x2_3, x1, x2_4, x1, x2_5, x1, x2_cf);
% obj(3).LineStyle = '--';
% obj(4).LineStyle = '--';
% obj(5).LineStyle = '--';
% obj(6).LineStyle = '--';
% obj(7).LineStyle = '--';
% obj(8).LineStyle = '--';
% % plot optimal solution:
% plot(valbis2(1), valbis2(2), 'r*');
% hold off % releases graph
% % set axes and figure legend:
% legend(obj, {'x \geq 0', 'y \geq 0', ...
%     'c1: y \geq 400 - x', ...
%     'c2: y \geq 600/2 - 1.3/2 x', ...
%     'c3: y \leq 1000/2 - 1.3/2 x', ...
%     'c4: y \geq 0.15/0.3x', ...
%     'c5: y \leq 3x', ...
%     'cf: y = -0.3x / 0.9'}, ...
%     'Location', 'Best');
% 
% [valbis2, fvalbis2] = linprog(f, Abis2, bbis2, [], [], lb, []);
% disp('val = ');
% disp(valbis2);
% disp('fval = ');
% disp(fvalbis2);
% 
% %% Sensitivity Analysis 3
% 
% % Cost function to minimize
% f = [0.30 0.90]; 
% 
% % New constraints matrix
% Abis3 = [-1     -1    ; 
%         -0.7   -2  ;
%          0.7    2  ;
%          0.15  -0.3  ;
%         -0.03   0.01 ;] ;
% 
% bbis3 = [-400 -600 1000 0 0]; 
% 
% figure(4)
% lb = [0 0]; % lower bounds
% ub = [1000 1000]; % upper bound (for the plot)
% [valbis3, fvalbis3] = linprog(f, Abis3, bbis3, [], [], lb, []);
% plotregion(-Abis3,-bbis3,lb,[],'g',0.5) % other file
% xlabel('x1: Corn'), ylabel('x2: Fodder')
% axis([lb(1) ub(1) lb(2) ub(2)]), grid
% title('Linear Problem for livestock feeding: Graph Solving')
% hold on
% 
% % Constraints lines:
% x1 = lb(1):1000; 
% x2 = lb(2):1000;
% x2_1 = 400 - x1; 
% x2_2 = 600/2 - 0.7/2*x1 ; 
% x2_3 = 1000/2 - 0.7/2*x1;
% x2_4 = 0.15/0.3*x1;
% x2_5 = 3*x1;
% 
% % Optimal cost function line:
% x2_cf = -(0.3.*x1) / 0.9; 
% obj = plot(x1, zeros(size(x1)), zeros(size(x2)), x2, x1, x2_1, x1, x2_2, x1, x2_3, x1, x2_4, x1, x2_5, x1, x2_cf);
% obj(3).LineStyle = '--';
% obj(4).LineStyle = '--';
% obj(5).LineStyle = '--';
% obj(6).LineStyle = '--';
% obj(7).LineStyle = '--';
% obj(8).LineStyle = '--';
% % plot optimal solution:
% plot(valbis3(1), valbis3(2), 'r*');
% hold off % releases graph
% % set axes and figure legend:
% legend(obj, {'x \geq 0', 'y \geq 0', ...
%     'c1: y \geq 400 - x', ...
%     'c2: y \geq 600/2 - 0.7/2 x', ...
%     'c3: y \leq 1000/2 - 0.7/2 x', ...
%     'c4: y \geq 0.15/0.3x', ...
%     'c5: y \leq 3x', ...
%     'cf: y = -0.3x / 0.9'}, ...
%     'Location', 'Best');
% 
% [valbis3, fvalbis3] = linprog(f, Abis3, bbis3, [], [], lb, []);
% disp('val = ');
% disp(valbis3);
% disp('fval = ');
% disp(fvalbis3);



% %% PART 3 
% 
% %% Problem
% 
% % Cost function to minimize
% f1 = [0.30 0.90 0.15 0.2]; 
% 
% % Constraint matrix
% A1 = [-1     -1    -2    -2; 
%      -1.3   -1.4  -2    -2;
%       1.3    1.4   2     2;
%       0.21  -0.3  -0.4  -0.6;
%      -0.03   0.01  0.1   0.08;
%       0      0     -1    -1] ;
% 
% b1 = [-400 -600 1000 0 0 2]; 
% 
% %% Solving by intlinprog
% 
% lb1 = [75, 75, 1, 1];
% ic = (3:4);
% [val1, fval1] = intlinprog(f1, ic, A1, b1, [], [], lb1, []);
% disp('val = ');
% disp(val1);
% disp('fval = ');
% disp(fval1);
% 
% 
% %% Sensitivity Analysis Integer 1
% 
% % Cost function to minimize
% f1 = [0.30 0.90 0.15 0.2]; 
% 
% % New constraints matrix
% A1bis = [-1     -1    -2    -2; 
%      -1.3   -1.4  -2    -2;
%       1.3    1.4   2     2;
%       0.15  -0.3  -0.4  -0.6;
%      -0.03   0.01  0.1   0.08;
%       0      0     -1    -1] ;
% 
% b1bis = [-400 -600 1000 0 0 2]; 
% 
% % Solving by intlinprog
% lb1 = [75, 75, 1, 1];
% ic = (3:4);
% [val1bis, fval1bis] = intlinprog(f1, ic, A1bis, b1bis, [], [], lb1, []);
% disp('val = ');
% disp(val1bis);
% disp('fval = ');
% disp(fval1bis);
% 
% %% Sensitivity Analysis Integer 2
% % Cost function to minimize
% f1 = [0.30 0.90 0.15 0.2]; 
%
% % New constraints matrix
% A1bis1 = [-1     -1    -2    -2; 
%           -1.3   -1.4  -3    -3;
%            1.3    1.4   3     3;
%            0.21  -0.3  -0.4  -0.6;
%           -0.03   0.01  0.1   0.08;
%            0      0     -1    -1] ;
% 
% b1bis1 = [-400 -600 1000 0 0 2]; 
% 
% % Solving by intlinprog
% lb1 = [75, 75, 1, 1];
% ic = (3:4);
% [val1bis1, fval1bis1] = intlinprog(f1, ic, A1bis1, b1bis1, [], [], lb1, []);
% disp('val = ');
% disp(val1bis1);
% disp('fval = ');
% disp(fval1bis1);
% 
% 
% %% Sensitivity Analysis Integer 3
% % Cost function to minimize
% f1 = [0.30 0.90 0.15 0.2]; 
%
% % New constraints matrix
% A1bis2 = [-1     -1    -2    -2; 
%      -1.3   -1.4  -2    -2;
%       1.3    1.4   2     2;
%       0.21  -0.3  -0.4  -0.6;
%      -0.03   0.01  0.1   0.08] ;
% 
% b1bis2 = [-400 -600 1000 0 0]; 
% 
% % Solving by intlinprog
% lb1bis2 = [75, 75, 0, 0];
% ic = (3:4);
% [val1bis2, fval1bis2] = intlinprog(f1, ic, A1bis2, b1bis2, [], [], lb1bis2, []);
% disp('val = ');
% disp(val1bis2);
% disp('fval = ');
% disp(fval1bis2);



%% PART 4

% Cost function to minimize
fun = @(x)0.3*x(1) + 0.9*x(2) + 0.15*x(3) + 0.2*x(4) + 0.5*x(5)^2 + 0.7*x(6)^2;

% Initial point
x0 = [75 75 1 1 1 1];

% Lower boundaries 
lb = [75 75 0 0 1 1];
ub = [400 400 50 50 5 5];

% x(3), x(4), x(5), x(6) are integers
intcon = [3 4 5 6];

nonlcon = @ellipsecons;
% ga_options = optimoptions('ga','PlotFcn', @gaplotbestf);
ga_options = optimoptions('ga');

% Running the Genetic Algorithm 40 times 
S = zeros(40,6);
F = zeros(1,40);
for i = 1:40
    [sol, fval] = ga(fun,6,[],[],[],[],lb,ub,nonlcon,intcon,ga_options);
    S(i,1:6) = [sol];
    F(1,i) = fval;
end

% Constraints 
function [c, ceq]= ellipsecons(x)
c(1) = -x(1)-x(2)-2*x(3)-2*x(4)-3*x(5)^2-3*x(6)^2+400;
c(2) = -1.3*x(1)-1.4*x(2)-2*x(3)-2*x(4)-6*x(5)^2-7*x(6)^2+600;
c(3) = 1.3*x(1)+1.4*x(2)+2*x(3)+2*x(4)+6*x(5)^2+7*x(6)^2-1000;
c(4) = 0.21*x(1)-0.3*x(2)-0.4*x(3)-0.6*x(4)+0.5*x(5)^2+0.4*x(6)^2;
c(5) = -0.03*x(1)+0.01*x(2)+0.1*x(3)+0.08*x(4)+0.08*x(5)^2+0.06*x(6)^2;
ceq = [];
end

% Sensitivity 1
% Cost function to minimize
fun = @(x)0.3*x(1)+0.9*x(2)+0.15*x(3)+0.2*x(4)+0.5*x(5)^2+0.7*x(6)^2;

% Initial point
x0 = [75 75 1 1 1 1];

% Lower boundaries 
lb = [75 75 0 0 1 1];
ub = [400 400 50 50 5 5];

% x(3), x(4), x(5), x(6) are integers
intcon = [3 4 5 6];

nonlcon = @ellipsecons;
ga_options = optimoptions('ga');%,'PlotFcn', @gaplotbestf);
% 'ga','PlotFcn', @gaplotbestf

% Running the GA 40 times 
S = zeros(40,6);
F = zeros(1,40);
for i = 1:40
    [sol, fval] = ga(fun,6,[],[],[],[],lb,ub,nonlcon,intcon,ga_options);
    S(i,1:6) = [sol];
    F(1,i) = fval;
end


% Constraints 
function [c, ceq]= ellipsecons(x)
c(1) = -x(1)-x(2)-2*x(3)-2*x(4)-4*x(5)^2-4*x(6)^2+400;
c(2) = -1.3*x(1)-1.4*x(2)-2*x(3)-2*x(4)-8*x(5)^2-9*x(6)^2+600;
c(3) = 1.3*x(1)+1.4*x(2)+2*x(3)+2*x(4)+8*x(5)^2+9*x(6)^2-1000;
c(4) = 0.21*x(1)-0.3*x(2)-0.4*x(3)-0.6*x(4)+0.5*x(5)^2+0.4*x(6)^2;
c(5) = -0.03*x(1)+0.01*x(2)+0.1*x(3)+0.08*x(4)+0.08*x(5)^2+0.06*x(6)^2;
ceq = [];
end

% Sensitivity 2
% Cost function to minimize
fun = @(x)0.3*x(1)+0.9*x(2)+0.15*x(3)+0.2*x(4)+0.5*x(5)^2+0.7*x(6)^2;

% Initial point
x0 = [75 75 1 1 1 1];

% Lower boundaries 
lb = [75 75 0 0 1 1];
ub = [400 400 50 50 5 5];

% x(3), x(4), x(5), x(6) are integers
intcon = [3 4 5 6];

nonlcon = @ellipsecons;
ga_options = optimoptions('ga');%,'PlotFcn', @gaplotbestf);
% 'ga','PlotFcn', @gaplotbestf

% Running the GA 40 times 
S = zeros(40,6);
F = zeros(1,40);
for i = 1:40
    [sol, fval] = ga(fun,6,[],[],[],[],lb,ub,nonlcon,intcon,ga_options);
    S(i,1:6) = [sol];
    F(1,i) = fval;
end

% Constraints 
function [c, ceq]= ellipsecons(x)
c(1) = -x(1)-x(2)-2*x(3)-2*x(4)-4*x(5)^2-4*x(6)^2+400;
c(2) = -1.3*x(1)-1.4*x(2)-2*x(3)-2*x(4)-8*x(5)^2-9*x(6)^2+600;
c(3) = 1.3*x(1)+1.4*x(2)+2*x(3)+2*x(4)+8*x(5)^2+9*x(6)^2-1000;
c(4) = 0.21*x(1)-0.3*x(2)-0.4*x(3)-0.6*x(4)+0.7*x(5)^2+0.6*x(6)^2;
c(5) = -0.03*x(1)+0.01*x(2)+0.1*x(3)+0.08*x(4)+0.12*x(5)^2+0.08*x(6)^2;
ceq = [];
end






