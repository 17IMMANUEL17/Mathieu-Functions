%  DESCRIPTION OF THE WORKSPACE CONTENTS
%  
% runelip          Primary interactive program showing
%                  modal animations and countour plots
%
% MembranePaper    Describes the analysis using Mathieu
%                  functions. To read the paper use the
%                  command open('MembranePaper.pdf')
%
%  Angular Mathieu Functions
%  ce,    se       ce(z,m,q), se(z,m,q)
%  cev,   sev      vectorized ce and se
%  cep,   sep      d/dz for ce and se
%  cepv,  sepv     vectorized ce and se
%
%  Radial Mathieu Functions
%  Mc1,   Ms1      Mc1(z,m,q), Ms1(z,m,q)    
%  Mc1v,  Ms1v     vectorized Mc1 and Ms1
%  Mc2,   Ms2      Mc2(z,m,q), Ms2(z,m,q)     
%  Mc2v,  Ms2v     vectorized Mc2 and Ms2
%  Mc1p,  Ms1p     d/dz of Mc1 and Ms1
%  Mc1pv, Ms1pv    vectorized Mc1p and Ms1p
%  Mc2p,  Ms2p     d/dz of Mc2 and Ms2
%  Mc2pv, Ms2pv    vectorized Mc2p and Ms2p
%  
%  Auxiliary Functions
%  bes             short name for besselj
%  beszeros        zeros of integer order
%                  Bessel functions
%  q2w, w2q        function to convert bet-
%                  q and w (frequency) in
%                  elliptical coordinates
%  matue           fundamental routine for
%                  eigenvalues and eigen-
%                  vectors defining the
%                  Mathieu functions
%  matuqrts        computes q-roots of 
%                  Mc1 and Ms1 to define
%                  elliptic membrane fre-
%                  quencies
%  plotmatu        plots ce,se,Mc1 and Ms1
%                  over a range of values
%  qplotmatu       plots Mathieu function
%                  dependence on q


%  Functions Producing Graphical Results
%  figure14_1 ...  these 14 functions generate
%  figure14_11     results matching Figures
%                  14.1 thru 14.11 of 'Compu-
%                  tation of Special Functions'
%                  by Zhang and Jin   
disp(' '), help listfunctions