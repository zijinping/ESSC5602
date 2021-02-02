% Prob 3.4

PC = 0.014;
PE = 1 - PC;

disp('probabilities in general population');
disp(sprintf('PC %f  PE %f  PC+PE %f', PC, PE, PC+PE));
disp('');

PYgC = 0.99;
PNgC = 1 - PYgC;
PYgE = 0.01;
PNgE = 1 - PYgE;

disp('conditional probabilities of test');
disp(sprintf('PYgC %f  PNgC %f  PYgC+PNgC %f', PYgC, PNgC, PYgC+PNgC));
disp(sprintf('PYgE %f  PNgE %f  PYgE+PNgE %f', PYgE, PNgE, PYgC+PNgC));
disp('');

PCgY = (PYgC * PC) / (PYgC*PC + PYgE * PE );

disp('Baysian probability');
disp(sprintf('PCgY %f', PCgY));
