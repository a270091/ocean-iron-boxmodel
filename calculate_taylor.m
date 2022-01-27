function [R,E,B,sf,sr] = calculate_taylor(refdata,moddata,weight)
%----------------------------------------------------------------
% function [R,E,B,sf,sr] = calculate_taylor(refdata,moddata,weight)
%
% calculates the correlation R, 
% the RMS-diference (after substracting the mean) E,
% the bias B (i.e. the difference of the means)
% and the standarddeviations SF, SR 
%
% between a 'reference data set' and a model output. The data can
% be weighted. Formulae are as given in 
%
% K.E. Taylor (2001): Summarizing multiple aspects of model
% performance in a single diagram. J. Geophys. Res. D, 106, 
% 7183-7192 
%----------------------------------------------------------------

  if (nargin<3),
    weight = ones(size(moddata));
  end
  
  % throw out pints where either model or reference have no data
  
  ii = find( ~isnan(refdata) & ~isnan(moddata) & ~isnan(weight));
  r = refdata(ii);
  f = moddata(ii);
  w = weight(ii);
  nd = length(ii);
  
  % sum of weight factors
  
  sw = sum(w(:));
  
  % weighted mean of data and model
  
  rm = sum( r(:).*w(:) ) / sw;
  fm = sum( f(:).*w(:) ) / sw;
  B = rm - fm;
  
  
  
  
  % variances of data and model
  
  rvar = sum( (r(:) - rm).^2 .*w(:) ) / sw;
  fvar = sum( (f(:) - fm).^2 .*w(:) ) / sw;
  
  % standard deviations
  
  sr = sqrt(rvar);
  sf = sqrt(fvar);
  
  % RMS diference after substraction of means
  
  E = sqrt( sum( ((r(:)-rm) - (f(:)-fm)).^2.*w(:) ) / sw );
  
  % correlation coefficient
  
  R = sum( ((r(:)-rm) .* (f(:)-fm)) .* w(:) ) / sw / sf / sr;
  
  return
  