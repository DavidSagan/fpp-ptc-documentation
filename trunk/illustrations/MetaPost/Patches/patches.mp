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
input geometrymarks
%input graph
%input gridticks
%input 3dgeom

%:my colors
%color yellow; yellow = red + green;
%color purple; purple = red + blue;
%color orange; orange = red + .5green;
color Hunter; Hunter = 0.3[DarkGreen,SteelBlue];
color gray;   gray = 0.95white;
% for ref: Crimson = (0.864, 0.080, 0.235);

%:local macros
numeric pt_size; pt_size := 3;
def drawpoint(expr z, c)=
  draw z withpen pencircle scaled pt_size withcolor c
enddef;

def drawsquare(expr s, o, c)=
  draw unitsquare scaled s shifted o withcolor c
enddef;

def drawrect(expr sx, sy, o, c)=
  draw unitsquare xscaled sx yscaled sy shifted o withcolor c
enddef;

def drawuvecs(expr s, o, ang, c)=
  drawarrow (origin -- s*dir(ang)) shifted o withcolor c;
  drawarrow (origin -- s*dir(ang+90)) shifted o withcolor c;
enddef;

secondarydef v projectedalong w =
  if pair(v) and pair(w):
    (v dotprod w) / (w dotprod w) * w
  else:
    errmessage("<projectedalong> :: arguments must be vectors")
  fi;
enddef;

string text_directive;
text_directive := "\footnotesize";
def mytext(expr s)=
  text_directive & "{" & s & "}"
enddef;


%:figure 1 :: patching
beginfig(1)
numeric u; u := 1.3cm;
numeric ang, el, ew, uvwt;
ang := -13;
el := 1.6u;
ew = 0.7; 
uvwt := 1.2;

picture elem;
elem := image(
  drawoptions(withpen pencircle scaled ew);
  drawrect(el,u, (0,-.5u), black);
  drawoptions(withpen pencircle scaled uvwt);
  drawuvecs(.35u, origin, 0, Hunter);
  drawuvecs(.35u, (el,0), 0, Hunter);
  drawoptions();
);

pair p[];
p0 := (el,0);
p1 := p0 + .5el*dir(ang);
p2 := p1 + el*dir(ang);
p3 := p0 + 2el*dir(ang);
p4 := (.5el,.7u);
p5 := p1 + (p4 rotated ang);
p6 := p3 + (p4 rotated ang);
p7 := .5[p5,p6];
p8 := .5[p0,p1] + .1u*dir(.5ang - 90);
p9 := p8 + .6u*dir(.5ang - 90);
p10 := .5[p2,p3] + .1u*dir(ang - 90);
p11 := p10 + .6u*dir(ang - 90);
p12 := (-.7el,0);
p13 := (-.1el,0);

draw p2 -- p3 withpen pencircle scaled .8 withcolor DarkCyan;
draw elem;
draw elem rotated ang shifted p1;
draw elem rotated ang shifted p3;

label(TEX(mytext("element 1")),p4);
label(TEX(mytext("element 2")),origin) rotated ang shifted p5;
label(TEX(mytext("element 4")),origin) rotated ang shifted p6;
label(TEX(mytext(
  "\begin{minipage}{1in}\centering "
  & "element 3\\ (drift)"
  & "\end{minipage}")),
  origin) rotated ang shifted p7;

drawarrow p12 -- p13;
label.urt(TEX(mytext(
  "\begin{minipage}{1cm} "
  & "particle\\ path"
  & "\end{minipage}")),
  origin) shifted p12;

drawoptions(withcolor DarkCyan);
draw p9 -- p8;
draw p11 -- p10;
label.bot(TEX(mytext(
  "\begin{minipage}{1in}\centering "
  & "patch\\ required"
  & "\end{minipage}")),
  origin) rotated .5ang shifted p9;
label.bot(TEX(mytext(
  "\begin{minipage}{1in}\centering "
  & "no patch\\ required"
  & "\end{minipage}")),
  origin) rotated ang shifted p11;
endfig;


%:figure 2 :: patching recirculator
beginfig(2)
numeric u; u := 1.3cm;
numeric ew, tw;
ew = 0.7; 
tw = 1.0;

numeric el, wd[];
el := 1.6u;
wd1 := 1.3u;
wd2 := 1u;

%% create bend
pair p[];
path bend;
-xpart p1 = xpart p2 =  .50el;
-xpart p4 = xpart p3 =  .75el;
 ypart p1 = ypart p2 = -.50wd1;
 ypart p3 = ypart p4 =  .50wd1;
bend := p1 -- p2 -- p3 -- p4 -- cycle;

%% create elem
picture elem;
elem := image(
  drawoptions();
  drawrect(el,wd2, (-.5el,-.5wd2), black);
);

%% create particle trajectories
numeric ang[], t[], d[];
d1 = 2.4el;
d2 = 2.0el;
pair c[], ctr[], e[];
path circ, traj[];
c0 = .5[p1,p4];
ctr0 = whatever[p1,p4] = whatever[p2,p3];
ctr1 = ctr0;
ctr2 = ctr0 - .8(p1 - p4);
for i = 1 upto 2:
  circ := fullcircle scaled 2arclength(ctr[i] -- c0) shifted ctr[i];
  (t[2i-1], whatever) = circ intersectiontimes (p1 -- p4);
  (t[2i], whatever) = circ intersectiontimes (p2 -- p3);
  ang[i] = angle(- direction t[2i] of circ);
  c[i] = point t[2i] of circ;
  e[i] = c[i] + d[i] * dir(ang[i]);
  traj[i] = (subpath (t[2i-1],t[2i]) of circ)
               -- (e[i] + (.5el+12pt)*dir(ang[i]));
endfor;

%% identify element shifts
c4 = e2 - .5el*dir(ang2);
c3 = whatever[c1,e1];
(c2 - c4) dotprod (c3 - c4) = 0;

%% draw marks
picture dist;
dist := dashpattern(on 0.7bp off 1.2bp);
c5 = c3 + .6u*dir(ang2 + 90);
c6 = c3 + .6u*dir(ang1 + 90);
c7 = c3 + .4u*dir(ang2);
drawoptions(withcolor DarkGoldenrod);
draw c4 -- c5;
draw c3 -- c6;
draw c3 -- c7;
mark_angle(c6, c3, c5, 2);
drawoptions(dashed dist);
drawdblarrow (c1 -- c3) shifted .6(c6 - c3) dashed dist;
drawdblarrow (c4 -- c3) shifted .5(c7 - c3) dashed dist;
label.urt(TEX(mytext("$d$")) rotated ang1, .5[c1,c3] shifted .6(c6 - c3));
label.rt(TEX(mytext("$h$")) rotated ang1, .65[c4,c3] shifted .6(c7 - c3));
label.rt(TEX(mytext("$\alpha$")) rotated ang1, .5[c3,c5]);

%% draw and label elements
drawoptions(withpen pencircle scaled ew);
draw bend;
for i = 1 upto 2:
  draw elem rotated ang[i] shifted e[i]; 
endfor;
label(TEX(mytext("magnet 1")), .35wd1*up);
label(TEX(mytext("magnet 2")), origin)
  rotated ang1 shifted (e1 + .3wd2*dir(ang1+90));
label(TEX(mytext("magnet 3")), origin)
  rotated ang2 shifted (e2 + .3wd2*dir(ang2+90));

%% draw trajectories
drawoptions(withpen pencircle scaled tw);
drawarrow traj1 withcolor Crimson;
drawarrow traj2 withcolor DarkCyan;

picture recirc;
recirc := currentpicture;
currentpicture := nullpicture;
draw recirc rotated -ang1;

endfig;


%:figure 3 :: misalignment
beginfig(3)
numeric u; u := 1.3cm;
numeric ang, dang, el, lite, uvwt;
ang := -13; dang := -6;
el := 1.6u;
lite := 0.6;
uvwt := 1.2;

picture elem, lightelem, redelem;
elem := image(
  drawrect(el,u, (0,-.5u), black);
  drawoptions(withpen pencircle scaled uvwt);
  drawuvecs(.35u, origin, 0, Hunter);
  drawuvecs(.35u, (el,0), 0, Hunter);
  drawoptions();
);

lightelem := image(
  drawrect(el,u, (0,-.5u), lite*white);
  drawoptions(withpen pencircle scaled uvwt);
  drawuvecs(.35u, origin, 0, lite[Hunter,white]);
  drawuvecs(.35u, (el,0), 0, lite[Hunter,white]);
  drawoptions();
);

redelem := image(
  drawrect(el,u, (0,-.5u), Crimson);
  drawoptions(withpen pencircle scaled uvwt);
  drawuvecs(.35u, origin, 0, Hunter);
  drawuvecs(.35u, (el,0), 0, Hunter);
  drawoptions();
);

pair p[];
p0 := (el,0);
p1 := p0 + .5el*dir(ang);
p2 := p1 + el*dir(ang);
p3 := p1 + ((.15u,.3u) rotated ang);
p4 := p3 + el*dir(ang+dang);
p5 := p0 + 2el*dir(ang);
p6 := (.5el,.7u);
p7 := p3 + (p6 rotated (ang + dang));
p8 := p5 + (p6 rotated ang);
p9 := (-.7el,0);
p10 := (-.1el,0);

draw elem;
draw lightelem rotated ang shifted p1;
draw redelem rotated (ang+dang) shifted p3;
draw elem rotated ang shifted p5;
drawarrow p1 -- p3 withcolor MediumBlue;
drawarrow p2 -- p4 withcolor MediumBlue;

label(TEX(mytext("element 1")),p6);
label(TEX(mytext("element 2")),origin)
  rotated (ang + dang) shifted p7 withcolor Crimson;
label(TEX(mytext("element 3")),origin)
  rotated ang shifted p8;

drawarrow p9 -- p10;
label.urt(TEX(mytext(
  "\begin{minipage}{1cm} "
  & "particle\\ path"
  & "\end{minipage}")),
  origin) shifted p9;

endfig;


end;
