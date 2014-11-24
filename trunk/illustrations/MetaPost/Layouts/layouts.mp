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
%input geometrymarks
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
color Hunter;     Hunter = 0.3[DarkGreen,SteelBlue];
color DarkYellow; DarkYellow = (1., 0.71, 0.16);
color gray;       gray = 0.95white;
% for ref: Crimson = (0.864, 0.080, 0.235);

%:local macros
numeric pt_size; pt_size := 3;
def drawpoint(expr z, c)=
  draw z withpen pencircle scaled pt_size withcolor c
enddef;

def square(expr s)=
  (-.5s,-.5s) -- (.5s,-.5s) -- (.5s,.5s) -- (-.5s,.5s) -- cycle
enddef;

def rect(expr sx, sy)=
  ((-.5sx,-.5sy) -- (.5sx,-.5sy) -- (.5sx,.5sy) -- (-.5sx,.5sy) -- cycle)
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
text_directive := "\footnotesize\textsf";
def mytext(expr s)=
  text_directive & "{" & s & "}"
enddef;


%:figure 1 :: simple ring
beginfig(1)
numeric u; u := 1.5cm;
numeric ew, pw;
ew = 0.5;
pw = 0.5;

numeric r, sx, sy;
r = .75u;
sx := 0.15u;
sy := 0.15u;

%% create magnet
picture elem;
elem := image(
  drawoptions(withpen pencircle scaled ew);
  unfill rect(sx,sy);
  draw rect(sx,sy);
  drawoptions();
);

%% draw ring
drawoptions(withpen pencircle scaled pw);
draw fullcircle scaled 2r;

%% place and label magnets
numeric a, da, nst;
nst = 8;
da  = 360 / nst;
for i = 0 upto 7:
  a := 90 - i * da;
  draw elem rotated a shifted (r*dir(a));
  label(textext(mytext("M" & decimal(i+1))), (r+10pt)*dir(a));
endfor;

endfig;


%:figure 2 :: recirculator
beginfig(2)
numeric u; u := 1.5cm;
numeric ew, pw;
ew = 0.5;
pw = 0.5;

numeric sx, sy;
sx := 0.15u;
sy := 0.15u;

%% create magnet
picture elem;
elem := image(
  drawoptions(withpen pencircle scaled ew);
  unfill rect(sx,sy);
  draw rect(sx,sy);
  drawoptions();
);
%% create cavity
picture cav;
cav := image(
  drawoptions(withpen pencircle scaled ew);
  unfill rect(10sx,sy);
  draw rect(10sx,sy);
  drawoptions();
);

%% create recirculator path
path recirc;
numeric c[], r[], dd;
-c1 = c2 = u;
dd = .68u;
r0 = 0.75u;
r1 = 0.5u;
r2 = 0.9u;
r3 = 1.3u;
r4 = 1.7u;
r5 = 2.1u;
pair p[];
ypart p1 = ypart p2  = ypart p3  = ypart p4  = ypart p5  = ypart p6 = r0;
ypart p7 = ypart p12 = ypart p13 = ypart p16 = ypart p17 = .6r0;
ypart p8 = ypart p11 = ypart p14 = ypart p15 = ypart p18 = -.6r0;
ypart p9 = ypart p10 = -r0;
xpart p1  = c1 - u - 4dd;
xpart p2  = c1 - u - 3dd;
xpart p3  = c1 - u - 2dd;
xpart p4  = c1 - u - 1dd;
xpart p5  = xpart p10 = c1;
xpart p6  = xpart p9  = c2;
xpart p7  = xpart p8  = c2 + r1;
xpart p12 = xpart p11 = c1 - r2;
xpart p13 = xpart p14 = c2 + r3;
xpart p16 = xpart p15 = c1 - r4;
xpart p17 = xpart p18 = c2 + r5;
recirc = p1
  for i = 2 upto 6: -- p[i]{right} endfor
  ... p[7]  .. p[8]  ...  {left}p[9] -- p[10]{left}
  ... p[11] .. p[12] ... {right}p[5] -- p[6]{right}
  ... p[13] .. p[14] ...  {left}p[9] -- p[10]{left}
  ... p[15] .. p[16] ... {right}p[5] -- p[6]{right}
  ... p[17] .. p[18] ...  {left}p[9] -- p[10]{left};

%% draw recirculator
drawoptions(withpen pencircle scaled pw);
draw recirc;
draw cav shifted ( r0*up);
draw cav shifted (-r0*up);

%% place and label magnets
numeric a, dl, di;
dl = 8pt;
pair pp;
%for i = 0 upto length(recirc):
%  a := angle(direction i of recirc);
%  draw elem rotated a shifted (point i of recirc);
%endfor;
for i = 1, 2, 3, 4, 5, 6:
  a := angle(direction i-1 of recirc);
  pp := point i-1 of recirc;
  draw elem rotated a shifted pp;
  label(textext(mytext("M" & decimal(i))), pp + dl*up);
endfor;
for i = 7, 8:
  a := angle(direction i-1 of recirc);
  pp := point i-1 of recirc;
  draw elem rotated a shifted pp;
  label(textext(mytext("M" & decimal(i))), pp + 1.5dl*left);
endfor;
for i = 9, 10:
  a := angle(direction i-1 of recirc);
  pp := point i-1 of recirc;
  draw elem rotated a shifted pp;
  label(textext(mytext("M" & decimal(i))), pp + dl*down);
endfor;
for i = 11, 12:
  a := angle(direction i-1 of recirc);
  pp := point i-1 of recirc;
  draw elem rotated a shifted pp;
  label(textext(mytext("M" & decimal(i))), pp + 1.6dl*right);
endfor;
for i = 13, 14:
  di := 2;
  a := angle(direction i-1+di of recirc);
  pp := point i-1+di of recirc;
  draw elem rotated a shifted pp;
  label(textext(mytext("M" & decimal(i))), pp + 1.6dl*left);
endfor;
for i = 15, 16:
  di := 4;
  a := angle(direction i-1+di of recirc);
  pp := point i-1+di of recirc;
  draw elem rotated a shifted pp;
  label(textext(mytext("M" & decimal(i))), pp + 1.6dl*right);
endfor;
for i = 17, 18:
  di := 6;
  a := angle(direction i-1+di of recirc);
  pp := point i-1+di of recirc;
  draw elem rotated a shifted pp;
  label(textext(mytext("M" & decimal(i))), pp + 1.6dl*left);
endfor;

endfig;


%:figure 3 :: recirculator, with fibres
beginfig(3)
numeric u; u := 1.5cm;
numeric ew, pw;
ew = 0.5;
pw = 0.5;

numeric sx, sy;
sx := 0.15u;
sy := 0.15u;

%% create magnet
picture elem;
elem := image(
  drawoptions(withpen pencircle scaled ew);
  unfill rect(sx,sy);
  draw rect(sx,sy);
  drawoptions();
);
%% create cavity
picture cav;
cav := image(
  drawoptions(withpen pencircle scaled ew);
  unfill rect(10sx,sy);
  draw rect(10sx,sy);
  drawoptions();
);
%% create fibre arrow
picture tinyarrow;
tinyarrow := image(
  ahlength := 3bp;
  ahangle := 45;
  drawoptions(withcolor Crimson);
  drawarrow 3left -- 3right;
);

%% create recirculator path
path recirc;
numeric c[], r[], dd;
-c1 = c2 = u;
dd = .68u;
r0 = 0.75u;
r1 = 0.5u;
r2 = 0.9u;
r3 = 1.3u;
r4 = 1.7u;
r5 = 2.1u;
pair p[];
ypart p1 = ypart p2  = ypart p3  = ypart p4  = ypart p5  = ypart p6 = r0;
ypart p7 = ypart p12 = ypart p13 = ypart p16 = ypart p17 = .6r0;
ypart p8 = ypart p11 = ypart p14 = ypart p15 = ypart p18 = -.6r0;
ypart p9 = ypart p10 = -r0;
xpart p1  = c1 - .9u - 4dd;
xpart p2  = c1 - .9u - 3dd;
xpart p3  = c1 - .9u - 2dd;
xpart p4  = c1 - .9u - 1dd;
xpart p5  = xpart p10 = c1;
xpart p6  = xpart p9  = c2;
xpart p7  = xpart p8  = c2 + r1;
xpart p12 = xpart p11 = c1 - r2;
xpart p13 = xpart p14 = c2 + r3;
xpart p16 = xpart p15 = c1 - r4;
xpart p17 = xpart p18 = c2 + r5;
recirc = p1
  for i = 2 upto 6: -- p[i]{right} endfor
  ... p[7]  .. p[8]  ...  {left}p[9] -- p[10]{left}
  ... p[11] .. p[12] ... {right}p[5] -- p[6]{right}
  ... p[13] .. p[14] ...  {left}p[9] -- p[10]{left}
  ... p[15] .. p[16] ... {right}p[5] -- p[6]{right}
  ... p[17] .. p[18] ...  {left}p[9] -- p[10]{left};

%% draw recirculator
drawoptions(withpen pencircle scaled pw);
draw recirc;
draw cav shifted ( r0*up);
draw cav shifted (-r0*up);

%% place and label magnets and fibres
numeric a, dl, di;
dl = 8pt;
pair pp, pf, pm, pa;
%% 1 .. 4
for i = 1 upto 4:
  a := angle(direction i-1 of recirc);
  pp := point i-1 of recirc;
  pm := pp + dl*up;
  pa := pm shifted 9left;
  pf := pa shifted left;
  draw elem rotated a shifted pp;
  text_directive := "\footnotesize\textsf";
  label(textext(mytext("M" & decimal(i))), pm);
  text_directive := "\scriptsize\textsf";
  label.lft(textext(mytext("F" & decimal(i))), pf) withcolor Crimson;
  draw tinyarrow shifted pa;
endfor;
%% 5
a := angle(direction 4 of recirc);
pp := point 4 of recirc;
pm := pp + dl*up;
pa := pm shifted 9left;
pf := pa shifted left;
draw elem rotated a shifted pp;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M5")), pm);
text_directive := "\scriptsize\textsf";
label.lft(textext(mytext("F5,\,F13,\,F21")), pf) withcolor Crimson;
draw tinyarrow shifted pa;
%% 6
a := angle(direction 5 of recirc);
pp := point 5 of recirc;
pm := pp + dl*up;
pa := pm shifted 9left;
pf := pa shifted left;
draw elem rotated a shifted pp;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M6")), pm);
text_directive := "\scriptsize\textsf";
label.lft(textext(mytext("F6,\,F14,\,F22")), pf) withcolor Crimson;
draw tinyarrow shifted pa;
%% 7 & 8
for i = 7, 8:
  a := angle(direction i-1 of recirc);
  pp := point i-1 of recirc;
  pm := pp + 1.5dl*left;
  pa := pm shifted 9left;
  pf := pa shifted left;
  draw elem rotated a shifted pp;
  text_directive := "\footnotesize\textsf";
  label(textext(mytext("M" & decimal(i))), pm);
  text_directive := "\scriptsize\textsf";
  label.lft(textext(mytext("F" & decimal(i))), pf) withcolor Crimson;
  draw tinyarrow shifted pa;
endfor;
%% 9
a := angle(direction 8 of recirc);
pp := point 8 of recirc;
pm := pp + dl*down;
pa := pm shifted 9left;
pf := pa shifted left;
draw elem rotated a shifted pp;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M9")), pm);
text_directive := "\scriptsize\textsf";
label.lft(textext(mytext("F9,\,F17,\,F25")), pf) withcolor Crimson;
draw tinyarrow shifted pa;
%% 10
a := angle(direction 9 of recirc);
pp := point 9 of recirc;
pm := pp + dl*down;
pa := pm shifted 11left;
pf := pa shifted left;
draw elem rotated a shifted pp;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M10")), pm);
text_directive := "\scriptsize\textsf";
label.lft(textext(mytext("F10,\,F18,\,F26")), pf) withcolor Crimson;
draw tinyarrow shifted pa;
%% 11 & 12
for i = 11, 12:
  a := angle(direction i-1 of recirc);
  pp := point i-1 of recirc;
  pm := pp + 1.6dl*right;
  pa := pm shifted 11right;
  pf := pa shifted right;
  draw elem rotated a shifted pp;
  text_directive := "\footnotesize\textsf";
  label(textext(mytext("M" & decimal(i))), pm);
  text_directive := "\scriptsize\textsf";
  label.rt(textext(mytext("F" & decimal(i))), pf) withcolor Crimson;
  draw tinyarrow rotated 180 shifted pa ;
endfor;
%% 13
a := angle(direction 14 of recirc);
pp := point 14 of recirc;
pm := pp + 1.6dl*left;
pa := pm shifted 7down;
pf := pa shifted down;
draw elem rotated a shifted pp;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M13")), pm);
text_directive := "\scriptsize\textsf";
label.bot(textext(mytext("F15")), pf) withcolor Crimson;
draw tinyarrow rotated 90 shifted pa;
%% 14
a := angle(direction 15 of recirc);
pp := point 15 of recirc;
pm := pp + 1.6dl*left;
pa := pm shifted 7up;
pf := pa shifted up;
draw elem rotated a shifted pp;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M14")), pm);
text_directive := "\scriptsize\textsf";
label.top(textext(mytext("F16")), pf) withcolor Crimson;
draw tinyarrow rotated -90 shifted pa;
%% 15
a := angle(direction 18 of recirc);
pp := point 18 of recirc;
pm := pp + 1.6dl*right;
pa := pm shifted 7up;
pf := pa shifted up;
draw elem rotated a shifted pp;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M15")), pm);
text_directive := "\scriptsize\textsf";
label.top(textext(mytext("F19")), pf) withcolor Crimson;
draw tinyarrow rotated -90 shifted pa;
%% 16
a := angle(direction 19 of recirc);
pp := point 19 of recirc;
pm := pp + 1.6dl*right;
pa := pm shifted 7down;
pf := pa shifted down;
draw elem rotated a shifted pp;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M16")), pm);
text_directive := "\scriptsize\textsf";
label. bot(textext(mytext("F20")), pf) withcolor Crimson;
draw tinyarrow rotated 90 shifted pa;
%% 17
a := angle(direction 22 of recirc);
pp := point 22 of recirc;
pm := pp + 1.6dl*left;
pa := pm shifted 7down;
pf := pa shifted down;
draw elem rotated a shifted pp;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M17")), pm);
text_directive := "\scriptsize\textsf";
label.bot(textext(mytext("F23")), pf) withcolor Crimson;
draw tinyarrow rotated 90 shifted pa;
%% 18
a := angle(direction 23 of recirc);
pp := point 23 of recirc;
pm := pp + 1.6dl*left;
pa := pm shifted 7up;
pf := pa shifted up;
draw elem rotated a shifted pp;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M18")), pm);
text_directive := "\scriptsize\textsf";
label.top(textext(mytext("F24")), pf) withcolor Crimson;
draw tinyarrow rotated -90 shifted pa;

text_directive := "\footnotesize\textsf";

endfig;


%% colors for the next two figures
color col[];
col1 := Black;
col2 := Teal;
col3 := Brown;
col4 := Green;
col5 := DarkOrchid;
col6 := .7[MediumBlue,RoyalBlue];
col7 := DarkGoldenrod;
col8 := VioletRed;


%:figure 4 :: DNA layouts
beginfig(4)
numeric u; u := 1.5cm;
numeric sx, sy;
sx := 0.38u;
sy := 0.20u;

pair l[], p[], m;
for i = 1 upto 8: l[i] := (0,-i*0.35u); endfor;
p0 := (.17u,0);
for i = 1 upto 4: p[i] := (i*0.54u,0); endfor;

picture tinyarrow;
tinyarrow := image(
  ahlength := 1.6bp;
  ahangle := 45;
  drawdblarrow (-.06u,0) -- (.06u,0);
);

drawoptions(withcolor col1);
label(textext(mytext("L1:")), l1 + p0);
m := l1 + p1;
%draw rect(sx,sy) shifted m;
draw tinyarrow shifted (m + .5(p2 - p1));
label(textext(mytext("M1")), m);
m := l1 + p2;
%draw rect(sx,sy) shifted m;
draw tinyarrow shifted (m + .5(p3 - p2));
label(textext(mytext("M2")), m);
m := l1 + p3;
%draw rect(sx,sy) shifted m;
draw tinyarrow shifted (m + .5(p4 - p3));
label(textext(mytext("M3")), m);
m := l1 + p4;
%draw rect(sx,sy) shifted m;
label(textext(mytext("M4")), m);

numeric im; im := 4;
string lbl;
for i = 2 upto 8:
  drawoptions(withcolor col[i]);
  lbl := mytext("L" & decimal(i) & ":");
  label(textext(lbl), l[i] + p0);
  im := im + 1;
  m := l[i] + p1;
  %draw rect(sx,sy) shifted m;
  draw tinyarrow shifted (m + .5(p2 - p1));
  lbl := mytext("M" & decimal(im));
  label(textext(lbl), m);
  im := im + 1;
  lbl := mytext("M" & decimal(im));
  m := l[i] + p2;
  %draw rect(sx,sy) shifted m;
  label(textext(lbl), m);
endfor;

endfig;


%:figure 5 :: recirculator, with fibres and layouts
beginfig(5)
numeric u; u := 1.5cm;
numeric ew, pw;
ew = 0.5;
pw = 0.5;

numeric sx, sy;
sx := 0.15u;
sy := 0.15u;

%% create magnet
%picture elem;
%color elemcolor;
%elemcolor := Black;
%elem := image(
%  drawoptions(withpen pencircle scaled ew);
%  unfill rect(sx,sy);
%  draw rect(sx,sy) withcolor elemcolor;
%  drawoptions();
%);
path elem;
elem := rect(sx,sy);
%% create cavity
%picture cav;
%cav := image(
%  drawoptions(withpen pencircle scaled ew);
%  unfill rect(10sx,sy);
%  draw rect(10sx,sy);
%  drawoptions();
%);
path cav;
cav := rect(10sx,sy);
%% create fibre arrow
picture tinyarrow;
tinyarrow := image(
  ahlength := 3bp;
  ahangle := 45;
  drawoptions(withcolor Crimson);
  drawarrow 3left -- 3right;
);

%% create recirculator path
path recirc;
numeric c[], r[], dd;
-c1 = c2 = u;
dd = .68u;
r0 = 0.75u;
r1 = 0.5u;
r2 = 0.9u;
r3 = 1.3u;
r4 = 1.7u;
r5 = 2.1u;
pair p[];
ypart p1 = ypart p2  = ypart p3  = ypart p4  = ypart p5  = ypart p6 = r0;
ypart p7 = ypart p12 = ypart p13 = ypart p16 = ypart p17 = .6r0;
ypart p8 = ypart p11 = ypart p14 = ypart p15 = ypart p18 = -.6r0;
ypart p9 = ypart p10 = -r0;
xpart p1  = c1 - .9u - 4dd;
xpart p2  = c1 - .9u - 3dd;
xpart p3  = c1 - .9u - 2dd;
xpart p4  = c1 - .9u - 1dd;
xpart p5  = xpart p10 = c1;
xpart p6  = xpart p9  = c2;
xpart p7  = xpart p8  = c2 + r1;
xpart p12 = xpart p11 = c1 - r2;
xpart p13 = xpart p14 = c2 + r3;
xpart p16 = xpart p15 = c1 - r4;
xpart p17 = xpart p18 = c2 + r5;
recirc = p1
  for i = 2 upto 6: -- p[i]{right} endfor
  ... p[7]  .. p[8]  ...  {left}p[9] -- p[10]{left}
  ... p[11] .. p[12] ... {right}p[5] -- p[6]{right}
  ... p[13] .. p[14] ...  {left}p[9] -- p[10]{left}
  ... p[15] .. p[16] ... {right}p[5] -- p[6]{right}
  ... p[17] .. p[18] ...  {left}p[9] -- p[10]{left};

%% draw recirculator
drawoptions(withpen pencircle scaled pw);
%draw recirc;
draw subpath ( 0, 4) of recirc withcolor col1;
draw subpath ( 4, 5) of recirc withcolor col2;
draw subpath ( 5, 8) of recirc withcolor col3;
draw subpath ( 8, 9) of recirc withcolor col4;
draw subpath ( 9,12) of recirc withcolor col5;
%draw subpath (12,13) of recirc withcolor col2;
draw subpath (13,16) of recirc withcolor col6;
%draw subpath (16,17) of recirc withcolor col4;
draw subpath (17,20) of recirc withcolor col7;
%draw subpath (20,21) of recirc withcolor col2;
draw subpath (21,24) of recirc withcolor col8;
%draw subpath (24,25) of recirc withcolor col4;
unfill cav shifted ( r0*up);
draw cav shifted ( r0*up) withcolor col2;
unfill cav shifted (-r0*up);
draw cav shifted (-r0*up) withcolor col4;

%% place and label magnets and fibres
numeric a, dl, di;
dl = 8pt;
pair pp, pf, pm, pa;
color elemcolor;
%% == layout 1 ==
elemcolor := col1;
%% 1 .. 4
for i = 1 upto 4:
  a := angle(direction i-1 of recirc);
  pp := point i-1 of recirc;
  pm := pp + dl*up;
  pa := pm shifted 9left;
  pf := pa shifted left;
  unfill elem rotated a shifted pp;
  draw elem rotated a shifted pp withcolor elemcolor;
  text_directive := "\footnotesize\textsf";
  label(textext(mytext("M" & decimal(i))), pm) withcolor elemcolor;
  text_directive := "\scriptsize\textsf";
  label.lft(textext(mytext("F" & decimal(i))), pf) withcolor Crimson;
  draw tinyarrow shifted pa withcolor Crimson;
endfor;
%% == layout 2 ==
elemcolor := col2;
%% 5
a := angle(direction 4 of recirc);
pp := point 4 of recirc;
pm := pp + dl*up;
pa := pm shifted 9left;
pf := pa shifted left;
unfill elem rotated a shifted pp;
draw elem rotated a shifted pp withcolor elemcolor;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M5")), pm) withcolor elemcolor;
text_directive := "\scriptsize\textsf";
label.lft(textext(mytext("F5,\,F13,\,F21")), pf) withcolor Crimson;
draw tinyarrow shifted pa withcolor Crimson;
%% 6
a := angle(direction 5 of recirc);
pp := point 5 of recirc;
pm := pp + dl*up;
pa := pm shifted 9left;
pf := pa shifted left;
unfill elem rotated a shifted pp;
draw elem rotated a shifted pp withcolor elemcolor;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M6")), pm) withcolor elemcolor;
text_directive := "\scriptsize\textsf";
label.lft(textext(mytext("F6,\,F14,\,F22")), pf) withcolor Crimson;
draw tinyarrow shifted pa withcolor Crimson;
%% == layout 3 ==
elemcolor := col3;
%% 7 & 8
for i = 7, 8:
  a := angle(direction i-1 of recirc);
  pp := point i-1 of recirc;
  pm := pp + 1.5dl*left;
  pa := pm shifted 9left;
  pf := pa shifted left;
  unfill elem rotated a shifted pp;
  draw elem rotated a shifted pp withcolor elemcolor;
  text_directive := "\footnotesize\textsf";
  label(textext(mytext("M" & decimal(i))), pm) withcolor elemcolor;
  text_directive := "\scriptsize\textsf";
  label.lft(textext(mytext("F" & decimal(i))), pf) withcolor Crimson;
  draw tinyarrow shifted pa withcolor Crimson;
endfor;
%% == layout 4 ==
elemcolor := col4;
%% 9
a := angle(direction 8 of recirc);
pp := point 8 of recirc;
pm := pp + dl*down;
pa := pm shifted 9left;
pf := pa shifted left;
unfill elem rotated a shifted pp;
draw elem rotated a shifted pp withcolor elemcolor;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M9")), pm) withcolor elemcolor;
text_directive := "\scriptsize\textsf";
label.lft(textext(mytext("F9,\,F17,\,F25")), pf) withcolor Crimson;
draw tinyarrow shifted pa withcolor Crimson;
%% 10
a := angle(direction 9 of recirc);
pp := point 9 of recirc;
pm := pp + dl*down;
pa := pm shifted 11left;
pf := pa shifted left;
unfill elem rotated a shifted pp;
draw elem rotated a shifted pp withcolor elemcolor;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M10")), pm) withcolor elemcolor;
text_directive := "\scriptsize\textsf";
label.lft(textext(mytext("F10,\,F18,\,F26")), pf) withcolor Crimson;
draw tinyarrow shifted pa withcolor Crimson;
%% == layout 5 ==
elemcolor := col5;
%% 11 & 12
for i = 11, 12:
  a := angle(direction i-1 of recirc);
  pp := point i-1 of recirc;
  pm := pp + 1.6dl*right;
  pa := pm shifted 11right;
  pf := pa shifted right;
  unfill elem rotated a shifted pp;
  draw elem rotated a shifted pp withcolor elemcolor;
  text_directive := "\footnotesize\textsf";
  label(textext(mytext("M" & decimal(i))), pm) withcolor elemcolor;
  text_directive := "\scriptsize\textsf";
  label.rt(textext(mytext("F" & decimal(i))), pf) withcolor Crimson;
  draw tinyarrow rotated 180 shifted pa withcolor Crimson;
endfor;
%% == layout 6 ==
elemcolor := col6;
%% 13
a := angle(direction 14 of recirc);
pp := point 14 of recirc;
pm := pp + 1.6dl*left;
pa := pm shifted 7down;
pf := pa shifted down;
unfill elem rotated a shifted pp;
draw elem rotated a shifted pp withcolor elemcolor;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M13")), pm) withcolor elemcolor;
text_directive := "\scriptsize\textsf";
label.bot(textext(mytext("F15")), pf) withcolor Crimson;
draw tinyarrow rotated 90 shifted pa withcolor Crimson;
%% 14
a := angle(direction 15 of recirc);
pp := point 15 of recirc;
pm := pp + 1.6dl*left;
pa := pm shifted 7up;
pf := pa shifted up;
unfill elem rotated a shifted pp;
draw elem rotated a shifted pp withcolor elemcolor;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M14")), pm) withcolor elemcolor;
text_directive := "\scriptsize\textsf";
label.top(textext(mytext("F16")), pf) withcolor Crimson;
draw tinyarrow rotated -90 shifted pa withcolor Crimson;
%% == layout 7 ==
elemcolor := col7;
%% 15
a := angle(direction 18 of recirc);
pp := point 18 of recirc;
pm := pp + 1.6dl*right;
pa := pm shifted 7up;
pf := pa shifted up;
unfill elem rotated a shifted pp;
draw elem rotated a shifted pp withcolor elemcolor;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M15")), pm) withcolor elemcolor;
text_directive := "\scriptsize\textsf";
label.top(textext(mytext("F19")), pf) withcolor Crimson;
draw tinyarrow rotated -90 shifted pa withcolor Crimson;
%% 16
a := angle(direction 19 of recirc);
pp := point 19 of recirc;
pm := pp + 1.6dl*right;
pa := pm shifted 7down;
pf := pa shifted down;
unfill elem rotated a shifted pp;
draw elem rotated a shifted pp withcolor elemcolor;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M16")), pm) withcolor elemcolor;
text_directive := "\scriptsize\textsf";
label. bot(textext(mytext("F20")), pf) withcolor Crimson;
draw tinyarrow rotated 90 shifted pa withcolor Crimson;
%% == layout 8 ==
elemcolor := col8;
%% 17
a := angle(direction 22 of recirc);
pp := point 22 of recirc;
pm := pp + 1.6dl*left;
pa := pm shifted 7down;
pf := pa shifted down;
unfill elem rotated a shifted pp;
draw elem rotated a shifted pp withcolor elemcolor;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M17")), pm) withcolor elemcolor;
text_directive := "\scriptsize\textsf";
label.bot(textext(mytext("F23")), pf) withcolor Crimson;
draw tinyarrow rotated 90 shifted pa withcolor Crimson;
%% 18
a := angle(direction 23 of recirc);
pp := point 23 of recirc;
pm := pp + 1.6dl*left;
pa := pm shifted 7up;
pf := pa shifted up;
unfill elem rotated a shifted pp;
draw elem rotated a shifted pp withcolor elemcolor;
text_directive := "\footnotesize\textsf";
label(textext(mytext("M18")), pm) withcolor elemcolor;
text_directive := "\scriptsize\textsf";
label.top(textext(mytext("F24")), pf) withcolor Crimson;
draw tinyarrow rotated -90 shifted pa withcolor Crimson;

text_directive := "\footnotesize\textsf";

endfig;


%:figure 6 :: RHIC main ring
beginfig(6)
numeric u; u := 1.3cm;

%% pen widths
numeric ewt, dwt;
ewt = 0.5;
dwt = 0.8;

%% RHIC
numeric dr[], r[];
dr1 = .5dwt;
dr2 = .04r0;
r0 = 2u;

%% set points on trajectories
pair ctr, ri[], ro[];
ri0 = r0 * right = ro0;
ri1 = ri0 + .30u * up + dr1 * left;
ri2 = ri0 + .45u * up + dr2 * left;
ro1 = ri1  + 2dr1 * right;
ro2 = ri2  + 2dr2 * right;
for i = 4 upto 6:
  ri[i] = ri[6-i] reflectedabout(origin, dir(30));
  ro[i] = ro[6-i] reflectedabout(origin, dir(30));
endfor;
ctr = ri2 + whatever*dir(0);
ctr = ri5 + whatever*dir(60);
r1 = arclength(ctr -- ri2);
r2 = arclength(ctr -- ro2);
ri3 = ctr + r1 * dir(30);
ro3 = ctr + r2 * dir(30);

%% define inner and outer trajectories
path rhic[];
rhic1 = (ri0 -- ri1 -- ri2{dir(90)} .. ri3 .. {dir(150)}ri4 -- ri5 -- ri6)
        rotated 30;
rhic2 = (ro0 -- ro1 -- ro2{dir(90)} .. ro3 .. {dir(150)}ro4 -- ro5 -- ro6)
        rotated 30;

%% set ring colors
color clr[];
clr1 := Orange;
clr2 := RoyalBlue;

%% draw RHIC
numeric da;
for i = 1 upto 6:
  da := i * 60;
  drawoptions(withpen pencircle scaled dwt withcolor clr1);
  draw rhic1 rotated da;
  drawoptions(withpen pencircle scaled dwt withcolor clr2);
  draw rhic2 rotated da;
  clr3 := clr1;
  clr1 := clr2;
  clr2 := clr3;
endfor;

%% zoom in on IP at 2'o'clock
numeric sx, sy;
sx := 1.5;
sy := 3;
numeric ex, ey;
color elemcolor;
ex := .15u;
ey := .15u;
elemcolor := Black;
picture elem[];
elem1 := image(
  drawoptions(withpen pencircle scaled ewt);
  unfill rect(ex,ey);
  draw rect(ex,ey);
  drawoptions();
);
elem2 := image(
  drawoptions(withpen pencircle scaled ewt);
  unfill rect(ex,5ey);
  draw rect(ex,5ey);
  drawoptions();
);
path ipstr;
picture ipzoom;
ipstr = ((ri0 -- ri1 -- ri2) shifted -ri0) xscaled sx yscaled sy;
pair ip[];
ip0 = sy * ypart(ri1) * up;
ip1 = (sy * ypart(ri2) + 2.5ey) * up + dr2 * sx * left;
ip2 = ip1 + 2dr2 * sx * right;
pair txshft;
txshft = 3pt * right;
ipzoom = image(
  drawoptions(withpen pencircle scaled dwt);
  draw ipstr withcolor RoyalBlue;
  draw ipstr reflectedabout(origin,up) withcolor Orange;
  draw ipstr reflectedabout(origin,right) withcolor Orange;
  draw ipstr rotated180 withcolor RoyalBlue;
  drawoptions();
  draw elem1 shifted ( 1.0ip0);
  draw elem1 shifted ( 0.6ip0);
  draw elem1 shifted (-0.6ip0);
  draw elem1 shifted (-1.0ip0);
  draw elem2 shifted ( ip1);
  draw elem2 shifted ( ip2);
  draw elem2 shifted (-ip2);
  draw elem2 shifted (-ip1);
  label.rt(textext(mytext("M1")), 1.0ip0 + txshft);
  label.rt(textext(mytext("M2")), 0.6ip0 + txshft);
  label.rt(textext(mytext("M3")),-0.6ip0 + txshft);
  label.rt(textext(mytext("M4")),-1.0ip0 + txshft);
);
%% draw it
draw ipzoom shifted (2r0 * right);
drawoptions(withcolor Hunter);
draw fullcircle scaled (2.03ypart(ri2)) shifted (ri0 rotated 30);
draw bbox ipzoom shifted (2r0 * right);
drawarrow ((ulcorner ipzoom shifted (2r0 * right)) -- (1.25(ri0 rotated 30)))
  shifted 2left;

endfig;


end;
