function kjm_printfig(fname,ppsize)
% This function exports the current figure in a reasonable way, in both eps
% and png formats kjm 05/10
%
% "fname" - desired filename. include path if desired
% "size" - size of figure in cm





set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [ppsize]);
set(gcf, 'PaperPosition',[0 0 2*ppsize])

print(gcf,fname,'-depsc2','-r300','-painters')
print(gcf,fname,'-dpng','-r300','-painters')


