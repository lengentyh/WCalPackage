ylabel('Energy (eV)');
xlabel('$\vec{k}/2\pi$','Interpreter','latex','FontSize',15);
xticks(sympoint);
xticklabels(sympointlabel);
set(gca,'xtick',xticks);
set(gca,'xticklabel',xticklabels);
axis([0 max(kpt(:,2))*normconstA 3 7]);
ax = gca;
ax.FontSize   = 18;
ax.FontWeight = 'bold';
ax.XGrid='on';

set(gcf, 'Position',  [150, 150, 2000, 1600])


