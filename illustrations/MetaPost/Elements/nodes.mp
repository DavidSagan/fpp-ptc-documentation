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
color DeepRed;  DeepRed = .5[Crimson,Maroon];
color Hunter;   Hunter = 0.3[DarkGreen,SteelBlue];
color gray;     gray = 0.90white;
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

string text_directive;
text_directive := "\small";
def mytext(expr s)=
  text_directive & "{" & s & "}"
enddef;

%% 3d macros
let vector = color;
let Xp = redpart;
let Yp = greenpart;
let Zp = bluepart;

vector vo,vi,vj,vk;
vo = (0,0,0);
vi = (1,0,0);
vj = (0,1,0);
vk = (0,0,1);

def sclprod(expr Va, Vb)= % scalar product
  (Xp(Va) * Xp(Vb) + Yp(Va) * Yp(Vb) + Zp(Va) * Zp(Vb))
enddef; 

def vecprod(expr Va, Vb)= % vector product
  (Yp(Va) * Zp(Vb) - Zp(Va) * Yp(Vb),
   Zp(Va) * Xp(Vb) - Xp(Va) * Zp(Vb),
   Xp(Va) * Yp(Vb) - Yp(Va) * Xp(Vb))
enddef;

def maxv(expr V)=
  max(abs(Xp(V)), abs(Yp(V)), abs(Zp(V)))
enddef;

def norm(expr V)=
  (maxv(V) * sqrt(sclprod(V / maxv(V), V / maxv(V))))
enddef;

def normed(expr V)=
  (V / norm(V))
enddef;

vardef rotateabout(expr Va, pt, axis, ang)=
  save v; vector v[];
  v0 = normed(axis);
  v1 = Va - pt;  % Va relative to pt
  v2 = sclprod(v1, v0) * v0; % part of v1 parallel to axis
  v3 = v1 - v2;              % part of v1 orthogonal to axis
  v4 = vecprod(v0,v3);
  v5 = pt + v2 + cosd(ang) * v3 + sind(ang) * v4;
  v5
enddef;

vardef arcabout(expr Va, pt, axis, ang, nst)=
  save a, da, p, v; vector v[];
  numeric a, da;
  da := ang / nst;
  
  v0 = normed(axis);
  v1 = Va - pt;  % Va relative to pt
  v2 = sclprod(v1, v0) * v0; % part of v1 parallel to axis
  v3 = v1 - v2;              % part of v1 orthogonal to axis
  v4 = vecprod(v0,v3);
  v5 = pt + v2 + cosd(ang) * v3 + sind(ang) * v4;
  v5
enddef;

vector observer_, screen_, vecup_, vecrt_;
vardef setprojection(expr o, s, u)=
  save v; vector v;
  observer_ := o;
  screen_ := s;
  v = whatever[o,u];
  sclprod(v-s,o-s) = 0;
  vecup_ := normed(v-s);
  vecrt_ := normed(vecprod(vecup_,o-s));
enddef;
setprojection(10vi + 3vj + 2vk, vo, vk);

vardef projection(expr P)=
  save Q, xp, yp;
  vector Q;
  Q = whatever[observer_,P];
  sclprod(Q-screen_,observer_-screen_) = 0;
  xp := sclprod(Q-screen_,vecrt_);
  yp := sclprod(Q-screen_,vecup_);
  (xp, yp)
enddef;


%:figure 1 :: bend with integration nodes
beginfig(1)
numeric ew, pw; % pen widths
ew := 2.5;
pw := 0.7;

numeric r[], rw;
rw := 2.5cm;
r0 := 4.5cm;
r1 := r0 - .5rw;
r2 := r0 + .5rw;

numeric ang[], ac, at, da, nst;
nst := 6;
ac := 90; % central angle
at := 70; % bend angle
da := at / nst;
ang0 := ac - .5at;
for i = 1 upto nst: ang[i] := ang0 + i * da; endfor;

path circ[], ray[];
circ1 := fullcircle scaled 2r1;
circ2 := fullcircle scaled 2r2;
for i = 0 upto nst: ray[i] := origin -- (1.2r2*dir(ang[i])); endfor;

%% draw body nodes
path wedge;
wedge = buildcycle(ray0,circ2,ray[nst],circ1);
fill wedge withcolor gray;
pair p[];
drawoptions(withpen pencircle scaled pw);
for i = 1 upto nst-1:
  p1 := ray[i] intersectionpoint circ1;
  p2 := ray[i] intersectionpoint circ2;
  draw p1 -- p2;
endfor;
draw wedge;

%% label body nodes
string mystr;
pair s;
drawoptions(withcolor Black);
mystr := "\begin{minipage}{4cm}\footnotesize\centering "
 & "any number ($N$) of body\\ integration nodes (here 6)"
 & "\end{minipage}";
label(TEX(mystr), (r2 + 18pt) * dir(ac));

%% draw and label exit nodes
p1 := ray0 intersectionpoint circ1;
p2 := ray0 intersectionpoint circ2;
%% exit fringe
drawoptions(withpen pencircle scaled ew withcolor DeepRed);
s := 0.6ew * dir(ang0 - 90);
draw (p1 -- p2) shifted s;
p3 := 1.02p2 shifted s;
p4 := p3 + .4rw*dir(ang0 + 25);
draw p3 -- p4 withpen pencircle scaled .5;
mystr := "\begin{minipage}{2.1cm}\footnotesize\centering "
 & "exit fringe\\ integration node"
 & "\end{minipage}";
label.top(TEX(mystr), p4);
%% exit patch
drawoptions(withpen pencircle scaled ew withcolor DarkCyan);
s := 1.6ew * dir(ang0 - 90);
draw (p1 -- p2) shifted s;
p3 := 1.02p2 shifted s;
p4 := p3 + .3rw*dir(ang0 - 10);
draw p3 -- p4 withpen pencircle scaled .5;
mystr := "\begin{minipage}{2.1cm}\footnotesize\centering "
 & "exit patch\\ integration node"
 & "\end{minipage}";
label.rt(TEX(mystr), p4);

%% draw and label entrance nodes
p1 := ray[nst] intersectionpoint circ1;
p2 := ray[nst] intersectionpoint circ2;
%% entrance fringe
drawoptions(withpen pencircle scaled ew withcolor DeepRed);
s := 0.6ew * dir(ang[nst] + 90);
draw (p1 -- p2) shifted s;
p3 := 1.02p2 shifted s;
p4 := p3 + .4rw*dir(ang[nst] - 25);
draw p3 -- p4 withpen pencircle scaled .5;
mystr := "\begin{minipage}{2.1cm}\footnotesize\centering "
 & "entrance fringe\\ integration node"
 & "\end{minipage}";
label.top(TEX(mystr), p4);
%% entrance patch
drawoptions(withpen pencircle scaled ew withcolor DarkCyan);
s := 1.6ew * dir(ang[nst] + 90);
draw (p1 -- p2) shifted s;
p3 := 1.02p2 shifted s;
p4 := p3 + .3rw*dir(ang[nst] + 10);
draw p3 -- p4 withpen pencircle scaled .5;
mystr := "\begin{minipage}{2.1cm}\footnotesize\centering "
 & "entrance patch\\ integration node"
 & "\end{minipage}";
label.lft(TEX(mystr), p4);

%% draw and label particle path
drawoptions();
p1 := .5[p1,p2] + 12pt*dir(ang[nst]+90);
p2 := p1 + 1.5cm*dir(ang[nst]+90);
drawarrow p2 -- p1;
mystr := "\begin{minipage}{0.9cm}\footnotesize "
 & "particle\\ path"
 & "\end{minipage}";
label.rt(TEX(mystr), origin)
  rotated (ang[nst] - 90)
  shifted (.0[p2,p1] + 12pt*dir(ang[nst]));

endfig;


%:figure 2 :: integration nodes with SC kick
beginfig(2)
numeric ew, pw; % pen widths
ew := 1.0;
pw := 0.4;
pt_size := 2.5;

numeric r[], rw;
rw := 4.0cm;
r0 := 6.0cm;
r1 := r0 - .5rw;
r2 := r0 + .5rw;

numeric ang[], ac, at, da, nst;
nst := 3;
ac := 105; % central angle
at := 55; % bend angle
da := at / nst;
ang0 := ac + .5at;
for i = 1 upto nst: ang[i] := ang0 - i * da; endfor;
ang[4] := ang1 - .4da;

path circ[], ray[];
circ1 := fullcircle scaled 2r1;
circ2 := fullcircle scaled 2r2;
for i = 0 upto nst+1: ray[i] := origin -- (1.2r2 * dir(ang[i])); endfor;

%% create nodes
path wedge;
wedge = buildcycle(ray0,circ2,ray[nst],circ1);
fill wedge withcolor gray;
pair p[];
drawoptions(withpen pencircle scaled pw);
for i = 1 upto nst-1:
  p1 := ray[i] intersectionpoint circ1;
  p2 := ray[i] intersectionpoint circ2;
  draw p1 -- p2;
endfor;
p1 := ray4 intersectionpoint circ1;
p2 := ray4 intersectionpoint circ2;
%% draw them
draw p1 -- p2 withcolor Crimson dashed evenly scaled .6;
draw wedge;
%% label planes
p1 := p1 shifted (-6pt*dir(ang[4]));
label(TEX("\small $s_\tau$"), p1) withcolor Crimson;
for i = 0 upto nst:
  p1 := ray[i] intersectionpoint circ1;
  p1 := p1 shifted (-6pt*dir(ang[i]));
  label(TEX("\small $s_" & decimal(i+1) & "$"), p1) withcolor Crimson;
endfor;

%% create path w/ no sc
numeric a[], tk;
pair pk;
path pno;
a0 := -0.60;
a1 := -0.05;
a2 :=  0.07;
a3 := -0.15;
for i = 0 upto 3: p[i] := (r0 + a[i]*.5rw) * dir(ang[i]); endfor;
pno := p0 for i = 1 upto 3: .. p[i] endfor;
p4 := p3 + .15rw*unitvector(direction 3 of pno);
%% draw it
drawoptions(withpen pencircle scaled ew withcolor Teal);
draw subpath (0,1) of pno;
draw subpath (1,3) of pno withcolor .5[Teal,white];
drawarrow p3 -- p4 withcolor .5[Teal,white];
%% drift before kick
picture drifting;
drifting := dashpattern(on 0.7bp off 1.2bp);
drawoptions(withcolor Teal);
pk := ray4 intersectionpoint (p1 -- (p1 + rw * (direction 1 of pno)));
draw p1 -- pk dashed drifting;
drawarrow pk -- (pk + .15rw*unitvector(pk-p1)) withcolor .5[Teal,white];
%% add labels
label.lrt(TEX("\footnotesize no SC kick"), p4);
drawpoint(p0, Teal);
drawpoint(p1, Teal);
drawpoint(p2, .5[Teal,white]);
drawpoint(p3, .5[Teal,white]);
label.lft(TEX("\small $t_1$"), p0);
label.lft(TEX("\small $t_2$"), p1);

%% create path w/ sc
path psc;
a1 := -0.25;
a2 :=  0.18;
a3 :=  0.15;
for i = 1 upto 3: p[i] := (r0 + a[i]*.5rw) * dir(ang[i]); endfor;
psc := p1{unitvector(pk-p1)} for i = 2 upto 3: .. p[i] endfor;
p4 := p3 + .15rw*unitvector(direction 3 of psc);
%% draw it
drawoptions(withpen pencircle scaled ew withcolor Crimson);
draw psc;
drawarrow p3 -- p4;
%% drift after kick
drawoptions(withcolor Crimson);
draw p1 -- pk dashed drifting;
drawarrow pk -- (pk + .15rw*unitvector(pk-p1));
%% add labels
label.rt(TEX("\footnotesize after SC kick"), p4);
drawpoint(pk, Crimson);
for i = 1 upto 3: drawpoint(p[i], Crimson); endfor;
label.lft(TEX("\small $\tau$"), pk);
label.llft(TEX("\small $t_2'$"), p1);
label.urt(TEX("\small $t_3$"), p2);
label.urt(TEX("\small $t_4$"), p3);

endfig;


end;
