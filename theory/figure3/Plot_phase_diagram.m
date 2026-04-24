
clear 
branch1=load('Phase_diagram_homogeneous.mat').branch1;
branch2=load('Phase_diagram_homogeneous.mat').branch2;

figure
plot(branch1(:,1),branch1(:,2),'k','linewidth',2)
hold on
plot(branch2(:,1),branch2(:,2),'k','linewidth',2)
axis square
xlabel('g');
ylabel('\nu_{ext}/\nu_{\theta}');
set(gca,'FontSize',18)
xlim([3 4.5])
ylim([0.5 2])
 