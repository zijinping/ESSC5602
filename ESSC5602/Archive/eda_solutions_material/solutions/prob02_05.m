% eda02_05
% load 'rocks.txt' and assign columns to mnemonic variable names
% file contaisn: SiO2	TiO2	Al2O3	FeOT	MgO	CaO	Na2O	K2O

disp('1=SiO2 2=TiO2 3=Al2O3 4=FeO_Total 5=MgO 6=CaO 7=Na2O 8=K2O');
disp('click to proceed to next plot');
disp('  click outside of axes to exit');

D = load('rocks.txt');
sio2 = D(:,1);   % SiO2
tio2 = D(:,2);   % TiO2
als03 = D(:,3);  % Al203
feot = D(:,4);   % FeO-total
mgo = D(:,5);    % MgO
cao = D(:,6);    % CaO
na20 = D(:,7);   % Na2O
k20 = D(:,8);    % K2O
Ns = size(D);
N = Ns(1);
M = Ns(2);

figure(1);
clf;
myplot=1;
for i = [1:M]
    Nbins=50;
    [myhist, mybins] = hist(D(:,i),Nbins);
    subplot(2,M/2,myplot);
    set(gca,'LineWidth',2);
    plot(mybins,myhist,'k-','LineWidth',2);
    xlabel(sprintf('element %d',i));
    ylabel('counts');
    myplot=myplot+1;
end


    
