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

input latexmp;
setupLaTeXMP(
  class = "tufte-handout",
  packages = "amsmath"
);

%:my colors
%color yellow; yellow = red + green;
%color purple; purple = red + blue;
%color orange; orange = red + .5green;
color Hunter; Hunter = 0.3[DarkGreen,SteelBlue];
color DarkYellow; DarkYellow = (1., 0.71, 0.16);
color gray;   gray = 0.90white;
color MediumGray; MediumGray = 0.72white;
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
text_directive := "\small\textsf";
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


%:figure 1 :: 3d "drift" block
beginfig(1)
numeric u; u := 1cm;
numeric el, lnwt, uvwt;
el := 1.0;
lnwt := 0.7;
uvwt := 1.2;

setprojection(25vi + 6vj + 6vk, vo, vk);

numeric sx, sy, sz;
sx := 4.5;
sy := 1.4;
sz := 1.0;

%% block
vector v[];
v0 := vo;
v1 := sy * vj;
v2 := sz * vk;
v3 := v1 + v2;
v4 := -sx * vi;
for i = 1 upto 3:
  v[i+4] := v[i] + v4;
endfor;
%% frames
v8 := .5[v0,v3];
v9  := v8 - el * vj;
v10 := v8 + el * vk;
v11 := .5[v4,v7];
v12 := v11 - el * vj;
v13 := v11 + el * vk;

%% projection
pair p[];
for i = 0 upto 13:
  p[i] := u * projection(v[i]);
endfor;
p14 := (p6 -- p7) intersectionpoint (p11 -- p13);

%% draw it
drawoptions(withpen pencircle scaled uvwt withcolor .5[Hunter,white]);
drawarrow p11 -- p12 dashed dashpattern(on 2.3bp off 2.8bp);
draw p11 -- p14;
drawarrow p14 -- p13 withcolor Hunter;

drawoptions(withpen pencircle scaled lnwt);
draw p5 -- p7 -- p6;
drawoptions(withcolor MediumGray);
draw p5 -- p4 -- p6;
draw p0 -- p4;
drawoptions(withpen pencircle scaled lnwt);
draw p1 -- p5;
draw p2 -- p6;
draw p3 -- p7;
draw p0 -- p1 -- p3 -- p2 -- cycle;

drawoptions(withpen pencircle scaled uvwt withcolor Hunter);
drawarrow p8 -- p9 dashed evenly;
drawarrow p8 -- p10;

drawoptions();
vector vr;
pair pr;
vr := .5[v2,v6];
pr := u * projection(vr);
%label.ulft(TEX(mytext("D")), pr);
label.ulft(textext(mytext("D")), pr);

endfig;


%:figure 2 :: 3d "bend" block
beginfig(2)
numeric u; u := 1cm;
numeric el, lnwt, uvwt;
el := 1.0;
lnwt := 0.7;
uvwt := 1.2;

setprojection(24vi + 5vj + 8vk, vo, vk);

numeric ang, rho, sx, sy, sz;
ang := -45;
rho := 3.0;
sx := 5.0;
sy := 1.4;
sz := 1.0;
vector axis, ctr;
axis := vk;
ctr := rho * vj;

%% block
vector v[];
v0 := vo;
v1 := sy * vj;
v2 := sz * vk;
v3 := v1 + v2;
for i = 0 upto 3:
  v[i+4] := rotateabout(v[i], ctr, axis, ang);
endfor;
%% frames
v8 := .5[v0,v3];
v9  := v8 - el * vj;
v10 := v8 + el * vk;
for i = 8 upto 10:
  v[i+3] := rotateabout(v[i], ctr, axis, ang);
endfor;
%% projection
pair p[];
for i = 0 upto 13:
  p[i] := u * projection(v[i]);
endfor;
p14 := u * projection(ctr);
p15 := (p6 -- p7) intersectionpoint (p11 -- p13);

%% draw it
numeric a, da, nst;
nst := 10;
da := ang / nst;
drawoptions(withpen pencircle scaled uvwt withcolor .5[Hunter,white]);
drawarrow p11 -- p12 dashed dashpattern(on 1.5bp off 2bp);
draw p11 -- p15;
drawarrow p15 -- p13 withcolor Hunter;

drawoptions(withpen pencircle scaled lnwt);
draw p5 -- p7 -- p6;
drawoptions(withcolor MediumGray);
draw p5 -- p4 -- p6;
%draw p0 .. p4;
draw p0 for i=1 upto nst:
  -- u * projection(rotateabout(v0, ctr, axis, i * da))
endfor;
drawoptions(withpen pencircle scaled lnwt);
%draw p1 .. p5;
%draw p2 .. p6;
%draw p3 .. p7;
for j=1 upto 3:
  draw p[j] for i=1 upto nst:
    -- u * projection(rotateabout(v[j], ctr, axis, i * da))
  endfor;
endfor;
draw p0 -- p1 -- p3 -- p2 -- cycle;

drawoptions(withpen pencircle scaled uvwt withcolor Hunter);
drawarrow p8 -- p9 dashed evenly;
drawarrow p8 -- p10;

drawoptions();
draw p1 -- p14 dashed evenly;
draw p5 -- p14 dashed evenly;

drawoptions();
vector vr;
pair pr;
vr := rotateabout(v2, ctr, axis, .53ang);
pr := u * projection(vr);
%label.ulft(TEX(mytext("B")), pr);
label.ulft(textext(mytext("B")), pr);

endfig;


%:figure 22 :: 3d "bend" block with s-axis
beginfig(22)
numeric u; u := .97cm;
numeric el, lnwt, uvwt;
el := 1.0;
lnwt := 0.7;
uvwt := 1.2;

setprojection(24vi + 5vj + 8vk, vo, vk);

numeric ang, rho, sx, sy, sz;
ang := -45;
rho := 3.0;
sx := 5.0;
sy := 1.4;
sz := 1.0;
vector axis, ctr;
axis := vk;
ctr := rho * vj;

%% block
vector v[];
v0 := vo;
v1 := sy * vj;
v2 := sz * vk;
v3 := v1 + v2;
for i = 0 upto 3:
  v[i+4] := rotateabout(v[i], ctr, axis, ang);
endfor;
%% frames
v8 := .5[v0,v3];
v9  := v8 - el * vj;
v10 := v8 + el * vk;
for i = 8 upto 10:
  v[i+3] := rotateabout(v[i], ctr, axis, ang);
endfor;
%% projection
pair p[];
for i = 0 upto 13:
  p[i] := u * projection(v[i]);
endfor;
p14 := u * projection(ctr + .5*(v2 - v0));
v15 := v11 + vecprod(v12-v11,v13-v11);
p15 := u * projection(v15);
p16 := (p5 -- p7) intersectionpoint (p11 -- p15);
p17 := (p6 -- p7) intersectionpoint (p11 -- p13);
p18 := (p5 -- p7) intersectionpoint (p11 -- p14);

%% draw it
numeric a, da, nst;
nst := 10;
da := ang / nst;

drawoptions(withpen pencircle scaled uvwt withcolor .5[Hunter,white]);
drawarrow p11 -- p12 dashed dashpattern(on 1.5bp off 2bp);
drawarrow p11 -- p13;

drawoptions(withpen pencircle scaled lnwt withcolor MediumGray);
draw p5 -- p4 -- p6;
%draw p0 .. p4;
draw p0 for i=1 upto nst:
  -- u * projection(rotateabout(v0, ctr, axis, i * da))
endfor;

%draw p8 .. p11;
drawoptions(withpen pencircle scaled uvwt withcolor .5[Hunter,white]);
draw p8 for i=1 upto nst:
  -- u * projection(rotateabout(v8, ctr, axis, i * da))
endfor;
draw p11 -- p16;
drawarrow p16 -- p15 withcolor Hunter;
drawarrow p17 -- p13 withcolor Hunter;

drawoptions(withpen pencircle scaled lnwt);
%draw p1 .. p5;
%draw p2 .. p6;
%draw p3 .. p7;
for j=1 upto 3:
  draw p[j] for i=1 upto nst:
    -- u * projection(rotateabout(v[j], ctr, axis, i * da))
  endfor;
endfor;
draw p0 -- p1 -- p3 -- p2 -- cycle;
draw p5 -- p7 -- p6;

drawoptions(withpen pencircle scaled uvwt withcolor Hunter);
drawarrow p8 -- p9 dashed evenly;
drawarrow p8 -- p10;

drawoptions(withpen pencircle scaled .5lnwt withcolor Crimson);
draw p8 -- p14;
draw p11 -- p18 withcolor .6[Crimson,white];
draw p18 -- p14;

drawoptions();
%label.lft(TEX("$x$"), p9);
%label.top(TEX("$y$"), p10);
%label.rt(TEX("$s$"), p15);
%label.bot(TEX("$\rho_0$"), .68[p8,p14] + 1up) withcolor Crimson;
%%label.urt(TEX("$\rho_0=1/\kappa_0$"), .6[p11,p14]) withcolor Crimson;
label.lft(textext(mytext("$x$")), p9);
label.top(textext(mytext("$y$")), p10);
label.rt(textext(mytext("$s$")), p15);
label.bot(textext(mytext("$\rho_0$")), .68[p8,p14] + 1up) withcolor Crimson;

endfig;


%:figure 3 :: "drift" w/ complex trajectory
beginfig(3)
numeric r[], dr, e;
r0 := 3cm;
dr := 0.25r0;

numeric ang[];
ang0 :=  55;
ang2 :=  20;
ang3 := ang0 - ang2;

r1 := 0.75r0;
r2 := r1/cosd(ang2);

pair p[];
p5 := r1*dir(ang0);
p6 := r2*dir(ang3);
currentpicture := nullpicture;

path rect;
e := arclength(p5--p6);
r3 := .28r0;
rect := unitsquare xscaled r0 yscaled 2r3;
fill rect withcolor gray;
draw rect withpen pencircle scaled 1bp;
draw (0,r3) -- (-e,r3) dashed evenly;
draw (r0,r3) -- (r0+e,r3) dashed evenly;

p0 := (0.00r0,1.3r3);
p1 := (0.20r0,1.2r3);
p2 := (0.30r0,0.8r3);
p3 := (0.50r0,0.8r3);
p4 := (0.60r0,1.2r3);
p5 := (0.80r0,1.4r3);
p6 := (1.00r0,0.3r3);
draw p0{dir(40)} for i=1 upto 5: .. p[i] endfor
     .. p6{dir(-05)}
  withpen pencircle scaled 2bp withcolor Crimson;
%for i=0 upto 6: drawpoint(p[i],yellow); endfor;

endfig;


%:figure 4 :: "bend" w/ complex trajectory
beginfig(4)
numeric r[];
r0 := 3cm;

numeric ang[];
ang0 :=  55;
ang1 := 125;
ang2 :=  20;
ang3 := ang0 - ang2;
ang4 := ang1 + ang2;

r1 := 0.75r0;
r2 := r1/cosd(ang2);

path circ[], ray[];
circ0 := fullcircle scaled 2r0;
circ1 := fullcircle scaled  r0;
ray0 := origin -- (1.1r0*dir(ang0));
ray1 := origin -- (1.1r0*dir(ang1));

path wedge;
wedge = buildcycle(ray0,circ0,ray1,circ1);
fill wedge withcolor gray;
draw wedge withpen pencircle scaled 1bp;

pair p[];
p0 := 1.10r1*dir(ang1);
p1 := 1.00r1*dir(110);
p2 := 0.90r1*dir(106);
p3 := 0.80r1*dir( 95);
p4 := 1.20r1*dir( 70);
p5 := 1.15r1*dir(ang0);
draw p0{dir(ang1-100)} for i=1 upto 4: .. p[i] endfor
     .. p5{dir(ang0-60)}
  withpen pencircle scaled 2bp withcolor Crimson;
%for i=0 upto 5: drawpoint(p[i],yellow); endfor;

p5 := r1*dir(ang0);
p6 := r2*dir(ang3);
p7 := r1*dir(ang1);
p8 := r2*dir(ang4);
draw p5 -- p6 dashed evenly;
draw p7 -- p8 dashed evenly;

endfig;


%:figure 34 :: "drift" and "bend" w/ complex trajectory
beginfig(34)
numeric r[], dr, e;
r0 := 3cm;
dr := 0.25r0;

numeric ang[];
ang0 :=  55;
ang1 := 125;
ang2 :=  20;
ang3 := ang0 - ang2;
ang4 := ang1 + ang2;

r1 := 0.75r0;
r2 := r1/cosd(ang2);

path circ[], ray[];
circ0 := fullcircle scaled 2r0;
circ1 := fullcircle scaled  r0;
ray0 := origin -- (1.1r0*dir(ang0));
ray1 := origin -- (1.1r0*dir(ang1));

path wedge;
wedge = buildcycle(ray0,circ0,ray1,circ1);
fill wedge withcolor gray;
draw wedge withpen pencircle scaled 1bp;

pair p[];
p0 := 1.10r1*dir(ang1);
p1 := 1.00r1*dir(110);
p2 := 0.90r1*dir(106);
p3 := 0.80r1*dir( 95);
p4 := 1.20r1*dir( 70);
p5 := 1.15r1*dir(ang0);
draw p0{dir(ang1-100)} for i=1 upto 4: .. p[i] endfor
     .. p5{dir(ang0-60)}
  withpen pencircle scaled 2bp withcolor Crimson;
%for i=0 upto 5: drawpoint(p[i],yellow); endfor;

p5 := r1*dir(ang0);
p6 := r2*dir(ang3);
p7 := r1*dir(ang1);
p8 := r2*dir(ang4);
draw p5 -- p6 dashed evenly;
draw p7 -- p8 dashed evenly;

picture bend;
bend := currentpicture;
currentpicture := nullpicture;

path rect;
e := arclength(p5--p6);
r3 := .28r0;
rect := unitsquare xscaled r0 yscaled 2r3;
fill rect withcolor gray;
draw rect withpen pencircle scaled 1bp;
draw (0,r3) -- (-e,r3) dashed evenly;
draw (r0,r3) -- (r0+e,r3) dashed evenly;

p0 := (0.00r0,1.3r3);
p1 := (0.20r0,1.2r3);
p2 := (0.30r0,0.8r3);
p3 := (0.50r0,0.8r3);
p4 := (0.60r0,1.2r3);
p5 := (0.80r0,1.4r3);
p6 := (1.00r0,0.3r3);
draw p0{dir(40)} for i=1 upto 5: .. p[i] endfor
     .. p6{dir(-05)}
  withpen pencircle scaled 2bp withcolor Crimson;
%for i=0 upto 6: drawpoint(p[i],yellow); endfor;

draw bend shifted (2.2r0,-1.5r3);

endfig;


%:figure 5 :: 3d connecting blocks
beginfig(5)
numeric u; u := .62cm;
numeric el, lnwt, uvwt;
el := 1.0;
lnwt := 0.8;
uvwt := 1.2;

vector vr;
pair pr;

setprojection(24vi + 5vj + 8vk, vo, vk);

numeric ang, rho, sx, sy, sz;
ang := -45;
rho := 3.0;
sx := 2.5;
sy := 1.4;
sz := 1.0;
vector axis[], ctr[];

%% blocks
vector v[];
v0 := vo;
v1 := sy * vj;
v2 := sz * vk;
v3 := v1 + v2;
ctr1 := rho * vj;
axis1 := vk;
for i = 0 upto 3:
  v[i+4] := rotateabout(v[i], ctr1, axis1, ang);
endfor;
vr = sx * normed(vecprod(v6 - v4, v5 - v4));
for i = 4 upto 7:
  v[i+4] := v[i] + vr;
endfor;
ctr2 := v8 + rho * vk;
axis2 := v9 - v8;
for i = 8 upto 11:
  v[i+4] := rotateabout(v[i], ctr2, axis2, -ang);
endfor;
%% frames
v16 := .5[v0,v3];
v17 := v16 - el * vj;
v18 := v16 + el * vk;
%% projection
pair p[];
for i = 0 upto 18:
  p[i] := u * projection(v[i]);
endfor;

%% draw it
numeric a, da, nst;
nst := 10;
da := ang / nst;

drawoptions(withcolor MediumGray);
draw p4 -- p8;
draw p5 -- p4 -- p6;
draw p9 -- p8 -- p10;
draw p13 -- p12 -- p14;
draw p0 for i=1 upto nst:
  -- u * projection(rotateabout(v0, ctr1, axis1, i * da))
endfor;
draw p8 for i=1 upto nst:
  -- u * projection(rotateabout(v8, ctr2, axis2, -i * da))
endfor;

drawoptions(withpen pencircle scaled lnwt);
draw p13 -- p15 -- p14;
draw p5 -- p7 -- p6;
draw p9 -- p11 -- p10;
draw p0 -- p1 -- p3 -- p2 -- cycle;
for i = 5 upto 7: draw p[i] -- p[i+4]; endfor;
for j=1 upto 3:
  draw p[j] for i=1 upto nst:
    -- u * projection(rotateabout(v[j], ctr1, axis1, i * da))
  endfor;
endfor;
for j=9 upto 11:
  draw p[j] for i=1 upto nst:
    -- u * projection(rotateabout(v[j], ctr2, axis2, -i * da))
  endfor;
endfor;

drawoptions(withpen pencircle scaled uvwt withcolor Hunter);
drawarrow p16 -- p17 dashed evenly;
drawarrow p16 -- p18;

drawoptions();
vr := rotateabout(v1, ctr1, axis1, .5ang);
pr := u * projection(vr);
%label.lrt(TEX(mytext("$\text{\textsf{B}}_\text{\textsf{H}}$")), pr);
label.lrt(textext(mytext("$\text{\textsf{B}}_\text{\textsf{H}}$")), pr);


vr := .4[v5,v9];
pr := u * projection(vr);
%label.lrt(TEX(mytext("D")), pr);
label.lrt(textext(mytext("D")), pr);
vr := rotateabout(v9, ctr2, axis2, -.5ang);
pr := u * projection(vr);
%label.lrt(TEX(mytext("$\text{\textsf{B}}_\text{\textsf{V}}$")), pr);
label.lrt(textext(mytext("$\text{\textsf{B}}_\text{\textsf{V}}$")), pr);

endfig;


%:figure 6 :: 3d connections
beginfig(6)
numeric u; u := .62cm;
numeric el, lnwt, uvwt;
el := 1.0;
lnwt := 0.8;
uvwt := 1.2;

vector vr;
pair pr;

setprojection(24vi + 5vj + 8vk, vo, vk);

numeric ang, rho, sx, sy, sz;
ang := -45;
rho := 3.0;
sx := 2.5;
sy := 1.4;
sz := 1.0;
vector axis[], ctr[];

%% blocks
vector v[];
v0 := vo;
v1 := sy * vj;
v2 := sz * vk;
v3 := v1 + v2;
ctr1 := rho * vj;
axis1 := vk;
for i = 0 upto 3:
  v[i+4] := rotateabout(v[i], ctr1, axis1, ang);
endfor;
vr = sx * normed(vecprod(v6 - v4, v5 - v4));
for i = 4 upto 7:
  v[i+4] := v[i] + vr;
endfor;
ctr2 := v8 + rho * vk;
axis2 := v9 - v8;
for i = 8 upto 11:
  v[i+4] := rotateabout(v[i], ctr2, axis2, -ang);
endfor;
%% frames
v16 := .5[v0,v3];
v17 := v16 - el * vj;
v18 := v16 + el * vk;
v19 := .5[v4,v7];
v20 := v19 - el * normed(v5 - v4);
v21 := v19 + el * vk;
v22 := .5[v8,v11];
v23 := v22 - el * vk;
v24 := v22 - el * normed(v9 - v8);
%% shrink block
for i = 4 upto 7: v[i+21] := .85[v[i],v[i+4]]; endfor;
%% block shift
vr := .5sx * normed(v8 - v4);
for i = 8 upto 15: v[i] := v[i] + vr; endfor;
for i = 22 upto 24: v[i] := v[i] + vr; endfor;
ctr2 := ctr2 + vr;
%% projection
pair p[];
for i = 0 upto 28:
  p[i] := u * projection(v[i]);
endfor;

%% draw it
numeric a, da, nst;
nst := 10;
da := ang / nst;

drawoptions(withcolor MediumGray);
draw p4 -- p25;
draw p26 -- p25 -- p27;
draw p13 -- p12 -- p14;
draw p8 for i=1 upto nst:
  -- u * projection(rotateabout(v8, ctr2, axis2, -i * da))
endfor;

drawoptions(withpen pencircle scaled lnwt);
draw p13 -- p15 -- p14;
draw p8 -- p9 -- p11 -- p10 -- cycle;
draw p26 -- p28 -- p27;
draw p4 -- p5 -- p7 -- p6 -- cycle;
for i = 5 upto 7: draw p[i] -- p[i+21]; endfor;
for j=9 upto 11:
  draw p[j] for i=1 upto nst:
    -- u * projection(rotateabout(v[j], ctr2, axis2, -i * da))
  endfor;
endfor;

drawoptions(withpen pencircle scaled uvwt withcolor Hunter);
drawarrow p19 -- p20 dashed dashpattern(on 1.5bp off 2bp);
drawarrow p19 -- p21;
drawarrow p22 -- p23 dashed dashpattern(on 1.5bp off 2bp);
drawarrow p22 -- p24;

drawoptions();
drawdblarrow (u * projection(v27 + .2vk)) -- (u * projection(v10 + .2vk));
vr := .6[v27, v10] + .2vk;
pr := u * projection(vr);
%label.ulft(TEX(mytext("?")), pr);
label.ulft(textext(mytext("?")), pr);

endfig;


end;
