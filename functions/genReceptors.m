% Draws initial TLR9 receptor levels from Gamma(2.2, 1420), the fit to the TLR9 flow data (Figure 5C).
function Ms1=genReceptors(mm)
rng default
% Ms1=gamrnd(4.32,1035,[mm 1]);   %TNF Receptor (Measured)
% Ms1=gamrnd(4,1014,[mm 1]);   %TLR2 Receptor (Measured)
Ms1=gamrnd(2.2,1420,[mm 1]);   %TLR9 Receptor (Measured)
% Ms2=gamrnd(2.2,1420,[mm 1]);   %TLR9 Receptor (Measured)
% Ms2=gamrnd(2.5,520,[mm 1]); 
% mean(Ms2)
% Ms1=Ms2;
% Ms3=gamrnd(4.80,4404,[mm 1]);   %TLR7 Receptor (Measured)
end