%!TEX TS-program = MPtoPDF

% use TEX("...") to create strings typeset by LaTeX
input TEX;
TEXPRE("%&latex"
       & char(10)
       & "\documentclass{tufte-handout}"
       & "\RequirePackage{amsmath}"
       & "\begin{document}");
TEXPOST("\end{document}");
%%label.lrt(TEX("${\cal A}^2$"),z2);

% filename template
% Allowed escape sequences for filenametemplate
%   %% A percent sign
%   %j The current jobname
%   %?0-9?c The charcode value
%   %?0-9?y The current year
%   %?0-9?m The numeric month
%   %?0-9?d The day of the month
%   %?0-9?H The hour
%   %?0-9?M The minute
%filenametemplate "%j.%c";  % default
filenametemplate "%j-%2c.eps";

% the next line
prologues:=3;
% means that
%   (a) MetaPost output will be EPSF, and
%   (b) fonts used will be embedded in the file
% for figures included in a tex document, prologues:=2 may suffice.

%warning check
%warningcheck:=0;

%:input packages
input colordef
%input graph
%input gridticks
%input 3dgeom

%:my colors
%color yellow; yellow = red + green;
%color purple; purple = red + blue;
%color orange; orange = red + .5green;
color Hunter; Hunter = 0.2[MediumSeaGreen,SteelBlue];
color gray;   gray = 0.85white;
% for ref: Crimson = (0.864, 0.080, 0.235);

%:local macros
numeric pt_size; pt_size := 3;
def drawpoint(expr z, c)=
  draw z withpen pencircle scaled pt_size withcolor c
enddef;

def drawsquare(expr s, o, c)=
  draw unitsquare scaled s shifted o withcolor c
enddef;

secondarydef v projectedalong w =
  if pair(v) and pair(w):
    (v dotprod w) / (w dotprod w) * w
  else:
    errmessage("<projectedalong> :: arguments must be vectors")
  fi;
enddef;


%:figure 1 :: one-turn maps on ring
beginfig(1)
numeric u; u := 1.2cm;
draw fullcircle scaled 2u;

%% radii
numeric r[];
r1 := 0.8;
r2 := 0.9;
r3 := 1.2;
%% directions
pair p[];
p1 := u*dir(80);
p2 := u*dir(30);
p3 := u*dir(-20);
%% map locations
for i=1 upto 3: draw r1*p[i] -- r3*p[i]; endfor;
%% ellipses [...]
pt_size := 1.5;
for i=-30 step -13 until -60: drawpoint((r2*u*dir(i)), black); endfor;
%% arrow 1--2
drawarrow r2*p1{dir(angle(p1)-90)} .. r2*p2{dir(angle(p2)-90)};

%% labels
label.top(TEX("$M_1$"),r3*p1);
label.rt(TEX("$M_2$"),r3*p2);
label.rt(TEX("$M_n$"),r3*p3);
label.lft(TEX("$M_{12}$"),r2*(.75[p1,p2]));

endfig;


%:figure 2 :: geometry of the normal form
beginfig(2)
numeric u; u := 1.0cm;
numeric w; w := 2mm;
numeric ang[], scl[];
pair p[], s[];
transform t[];

%% wheel
picture wheel;
fill unitsquare scaled .5w;
fill unitsquare scaled .5w shifted -.5(w,w);
clip currentpicture to (fullcircle scaled w);
draw fullcircle scaled w withpen pencircle scaled .3;
wheel := currentpicture;

%% ellipse 1
ang1 := 25;
scl1 := 2;
p1 := dir(15);
fill fullcircle xscaled (scl1*u) yscaled (u/scl1) rotated ang1 withcolor gray;
draw fullcircle xscaled (scl1*u) yscaled (u/scl1) rotated ang1;
draw wheel shifted (p1 xscaled (scl1*u/2) yscaled (u/scl1/2) rotated ang1);
%% circle 1
s1 := (5u,0);
fill fullcircle scaled u shifted s1 withcolor gray;
draw fullcircle scaled u shifted s1;
draw wheel shifted (p1 scaled u/2 shifted s1);

%% ellipse 2
ang2 := 70;
scl2 := 1.6667;
p2 := dir(-75);
s2 := (0,-3u);
fill fullcircle xscaled (scl2*u) yscaled (u/scl2) rotated ang2 shifted s2 withcolor gray;
draw fullcircle xscaled (scl2*u) yscaled (u/scl2) rotated ang2 shifted s2;
draw wheel shifted (p2 xscaled (scl2*u/2) yscaled (u/scl2/2) rotated ang2 shifted s2);
%% circle 2
s3 := s1 + s2;
fill fullcircle scaled u shifted s3 withcolor gray;
draw fullcircle scaled u shifted s3;
draw (origin -- .7u*p1) shifted s3 withcolor Crimson;
draw (origin -- .7u*p2) shifted s3 withcolor Crimson;
drawarrow ((.65u*p1){dir(angle(p1)-90)} .. (.65u*p2){dir(angle(p2)-90)}) shifted s3
  withpen pencircle scaled .4 withcolor Crimson;
draw wheel shifted (p1 scaled u/2 shifted s3) withcolor .5white;
draw wheel shifted (p2 scaled u/2 shifted s3);

%% arrows
p3 := 0.5u*dir( -10);
p4 := 0.8u*dir(  55) shifted s2;
p5 := 1.1u*dir(  25);
p6 := 0.6u*dir( 135) shifted s1;
p7 := 0.6u*dir( -30) shifted s1;
p8 := 0.6u*dir(  30) shifted s3;
p9 := 0.6u*dir(-135) shifted s3;
p10:= 0.4u*dir( -40) shifted s2;
drawarrow p3{dir( -55)} .. p4 {dir(-125)};
drawarrow p5{dir(  25)} .. p6 {dir( -25)};
drawarrow p7{dir( -55)} .. p8 {dir(-125)};
drawarrow p9{dir(-155)} .. p10{dir( 155)};

%% labels
label.lft(TEX("$M_{12}$"), .5[p3,p4] shifted (.3u,0));
label.bot(TEX("$A_1^{-1}$"), .5[p5,p6] shifted (0,.3u));
label.lft(TEX("$R(\Delta\phi_{12})$"), .5[p7,p8] shifted (.3u,0));
label.top(TEX("$A_2$"), .5[p9,p10] shifted (0,-.3u));
label.top(TEX("$M_{12} = A_2 \circ R(\Delta\phi_{12}) \circ A_1^{-1}$"),
          .5[p5,p6] shifted (0,.5u));
label.rt(TEX("$\Delta\phi_{12}$"), .65u*dir(-30) shifted s3) withcolor Crimson;

endfig;


end;
