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
input colordef;
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
color DeepRed;  DeepRed = .5[Crimson,Maroon];
color Hunter;   Hunter = 0.3[DarkGreen,SteelBlue];
color gray;     gray = 0.95white;
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

vardef ngon(expr n)=
  save i; numeric i;
  ((1,0) for i = 1 upto n-1: -- dir(i*360/n) endfor -- cycle)
enddef;

def drawsquare(expr s, o, c)=
  draw unitsquare scaled s shifted o withcolor c
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


%:global definitions
numeric u; u := 1cm;

%% pen width
numeric pw, dwt, twt, uvwt;
pw := 0.4;
dwt := 0.6;
twt := 0.6;
uvwt := 1.2;
%% dashing
picture pdots;
pdots = dashpattern (on 1pt off 1.5pt);

%% accel geometry
numeric da, ncell, rho;
ncell := 10;
da = 360 / ncell;
rho := 1.9u;
numeric ew, lb, lq;
ew := .25u;
lb := .40u;
lq := .08u;
pair ctr[], e[];
numeric rb;
e2 = -e1 = .5(lb,0);
ctr0 = e1 + whatever*dir(.5da - 90);
ctr0 = e2 + whatever*dir(-.5da - 90);
rb = arclength(ctr0 -- e1);
numeric esag, sag;
sag = 2rb * (sind(.25da) ** 2);
esag = sag / cosd(.5da);
numeric aq, rm, rp, rq;
aq := 9.5;
rp = rho + .5sag + esag;
rm = rp * cosd(.5da) = rq * cosd(.5da - aq);
ctr1 = ctr0 + .5sag*down;
path circ[];
numeric t[];
circ0 = fullcircle scaled 2rb shifted ctr1;
(t1, whatever) = circ0 intersectiontimes (.5(-lb,-ew) -- .5(-lb,ew));
(t2, whatever) = circ0 intersectiontimes (.5( lb,-ew) -- .5( lb,ew));
path sb;
numeric dl;
%dl := ew * sind(.5da) / cosd(.5da);
%sb = (-.5lb,-.5ew) -- (.5lb,-.5ew) -- (.5lb+dl,.5ew) -- (-.5lb-dl,.5ew)
%     -- cycle;
path ray[];
circ1 = fullcircle scaled (2rb + ew) shifted ctr1;
circ2 = fullcircle scaled (2rb - ew) shifted ctr1;
ray1 = (origin -- 2rho * dir(90 - .5da)) shifted ctr1;
ray2 = (origin -- 2rho * dir(90 + .5da)) shifted ctr1;
sb = buildcycle(ray1, circ1, ray2, circ2);


%% create elements
picture arcp, bend, quad, sbend, strb, tic;
arcp := image(
  drawoptions(withpen pencircle scaled pw withcolor Crimson);
  draw subpath (t1,t2) of circ0 dashed pdots;
  drawoptions();
);
bend := image(
  drawoptions(withpen pencircle scaled pw withcolor Crimson);
  fill rect(lb,ew) withcolor gray;
  draw rect(lb,ew);
  draw subpath (t1,t2) of circ0 dashed pdots;
  drawoptions();
);
quad := image(
  drawoptions(withpen pencircle scaled pw withcolor DarkBlue);
  fill rect(lq,ew) withcolor gray;
  draw rect(lq,ew);
  drawoptions();
);
sbend := image(
  drawoptions(withpen pencircle scaled pw withcolor Crimson);
  fill sb withcolor gray;
  draw sb;
  draw subpath (t1,t2) of circ0 dashed pdots;
  drawoptions();
);
strb := image(
  drawoptions(withpen pencircle scaled pw withcolor Hunter);
  fill rect(lb,ew) withcolor gray;
  draw rect(lb,ew);
  draw (-.5lb,0) -- (.5lb,0) dashed pdots;
  drawoptions();
);
tic := image(
  drawoptions(withpen pencircle scaled pw withcolor DarkCyan);
  draw (-.2ew,0) -- (.2ew,0);
  drawoptions();
);

%% create direction arrows
picture arrcw, arrccw;
arrcw := image(
  drawoptions(withpen pencircle scaled pw withcolor DarkCyan);
  drawarrow (.75rm * dir(90 + .5da)){dir(.5da)}
            .. {dir(-.5da)}(.75rm * dir(90 - .5da));
  drawoptions();
);
arrccw = arrcw xscaled -1;

%% create layouts
numeric ang;
picture psr, psr_minus, psr_qfb, psr_sb, psr_sb_minus;
psr = image(
  drawoptions(withpen pencircle scaled dwt);
  draw ngon(10) scaled rp;
  for i = 0 upto ncell - 1:
    ang := i * da;
    drawoptions();
    draw bend rotated (ang - 90) shifted (rho*dir(ang));
    draw quad rotated (ang - 90 + .5da) shifted (rq*dir(ang + aq));
    draw quad rotated (ang - 90 - .5da) shifted (rq*dir(ang - aq));
    drawoptions(withpen pencircle scaled twt);
    draw tic rotated (ang + .5da) shifted (rm*dir(ang + .5da));
  endfor;
  drawoptions();
);
psr_minus = image(
  drawoptions(withpen pencircle scaled dwt);
  %draw ngon(10) scaled rp;
  draw ((1,0) for i = 1 upto 9: -- dir(i*360/10) endfor)
    rotated -2da scaled rp;
  for i = -1 upto 6:
    ang := i * da;
    drawoptions();
    draw bend rotated (ang - 90) shifted (rho*dir(ang));
    draw quad rotated (ang - 90 + .5da) shifted (rq*dir(ang + aq));
    draw quad rotated (ang - 90 - .5da) shifted (rq*dir(ang - aq));
    drawoptions(withpen pencircle scaled twt);
    draw tic rotated (ang + .5da) shifted (rm*dir(ang + .5da));
  endfor;
  drawoptions();
  draw bend rotated (-2da - 90) shifted (rho*dir(-2da));
  draw quad rotated (-2da - 90 + .5da) shifted (rq*dir(-2da + aq));
  draw bend rotated (7da - 90) shifted (rho*dir(7da));
  draw quad rotated (7da - 90 - .5da) shifted (rq*dir(7da - aq));
  drawoptions(withpen pencircle scaled twt);
  draw tic rotated (-2da + .5da) shifted (rm*dir(-2da + .5da));
  drawoptions();
);
psr_qfb = image(
  for i = 0 upto ncell - 1:
    ang := i * da;
    draw strb rotated (ang - 90) shifted (rho*dir(ang));
  endfor;
  drawoptions();
);
psr_sb = image(
  drawoptions(withpen pencircle scaled dwt);
  draw ngon(10) scaled rp;
  for i = 0 upto ncell - 1:
    ang := i * da;
    drawoptions();
    draw sbend rotated (ang - 90) shifted (rho*dir(ang));
    draw quad rotated (ang - 90 + .5da) shifted (rq*dir(ang + aq));
    draw quad rotated (ang - 90 - .5da) shifted (rq*dir(ang - aq));
    drawoptions(withpen pencircle scaled twt);
    draw tic rotated (ang + .5da) shifted (rm*dir(ang + .5da));
  endfor;
  drawoptions();
);
psr_sb_minus = image(
  drawoptions(withpen pencircle scaled dwt);
  %draw ngon(10) scaled rp;
  draw ((1,0) for i = 1 upto 9: -- dir(i*360/10) endfor)
    rotated -2da scaled rp;
  for i = -1 upto 6:
    ang := i * da;
    drawoptions();
    draw sbend rotated (ang - 90) shifted (rho*dir(ang));
    draw quad rotated (ang - 90 + .5da) shifted (rq*dir(ang + aq));
    draw quad rotated (ang - 90 - .5da) shifted (rq*dir(ang - aq));
    drawoptions(withpen pencircle scaled twt);
    draw tic rotated (ang + .5da) shifted (rm*dir(ang + .5da));
  endfor;
  drawoptions();
  draw sbend rotated (-2da - 90) shifted (rho*dir(-2da));
  draw quad rotated (-2da - 90 + .5da) shifted (rq*dir(-2da + aq));
  draw tic rotated (-2da + .5da) shifted (rm*dir(-2da + .5da));
  draw sbend rotated (7da - 90) shifted (rho*dir(7da));
  draw quad rotated (7da - 90 - .5da) shifted (rq*dir(7da - aq));
  drawoptions();
);

%% points for placing layouts
pair p[];
p1 = 2(rp + ew) * left;
p2 = 2(rp + ew) * right;
p3 = (rm + 14pt) * down;
p4 = 2rm * up;
p5 = rm * up;


%:figure 1: siamese linked-list
beginfig(1)
numeric sl, sw;
sl := 13mm;
sw :=  4mm;

%% siamese boxes
pair sh;
sh = 1.8sw*up;
drawoptions(withpen pencircle scaled dwt);
draw rect(sl, sw);
draw rect(sl, sw) shifted sh;
drawoptions();

%% arrows
interim ahlength := 9bp;
interim ahangle := 15;

pair a[];
a1 = .28sl*left;
a2 = -a1;
a3 = a1 + sh;
a4 = a2 + sh;
drawarrow a1 -- a3;
drawarrow a4 -- a2;

%% labels
text_directive := "\footnotesize";
dotlabel.lft(textext(mytext("s")), a1);
dotlabel.rt(textext(mytext("s")), a4);

endfig;


%:figure 2: siamese linked-list, rotated incorrectly
beginfig(2)
numeric sl, sw, ang;
sl := 13mm;
sw :=  4mm;
ang := -20;

%% siamese boxes
pair sh;
sh = 1.8sw*up;
drawoptions(withpen pencircle scaled dwt);
draw rect(sl, sw) rotated ang;
draw rect(sl, sw) rotated ang shifted sh;
drawoptions();

%% arrows
interim ahlength := 9bp;
interim ahangle := 15;

pair a[];
a1 = .28sl*left rotated ang;
a2 = -a1;
a3 = a1 + sh;
a4 = a2 + sh;
drawarrow a1 -- a3;
drawarrow a4 -- a2;

%% labels
text_directive := "\footnotesize";
dotlabel.lft(textext(mytext("s")), a1);
dotlabel.rt(textext(mytext("s")), a4);

endfig;


%:figure 3: siamese linked-list, rotated correctly
beginfig(3)
numeric sl, sw, ang;
sl := 13mm;
sw :=  4mm;
ang := -20;

%% siamese boxes
pair sh;
sh = 1.8sw*up;
drawoptions(withpen pencircle scaled dwt);
draw rect(sl, sw) rotated ang;
draw rect(sl, sw) shifted sh rotated ang;
drawoptions();

%% arrows
interim ahlength := 9bp;
interim ahangle := 15;

pair a[];
a1 = .28sl*left rotated ang;
a2 = -a1;
a3 = a1 + sh rotated ang;
a4 = a2 + sh rotated ang;
drawarrow a1 -- a3;
drawarrow a4 -- a2;

%% labels
text_directive := "\footnotesize";
dotlabel.lft(textext(mytext("s")), a1);
dotlabel.rt(textext(mytext("s")), a4);

endfig;


%:figure 4: girder linked-list
beginfig(4)
numeric gl, gw;
gl := 8mm;
gw :=  4mm;

%% girder boxes
pair sh;
sh = 1.5gl*right;
drawoptions(withpen pencircle scaled dwt);
draw rect(gl, gw);
draw rect(gl, gw) shifted 1sh;
draw rect(gl, gw) shifted 2sh;
drawoptions();

%% points
pair a[];
a1 = .3gl*right;
a2 = a1 + 1sh;
a3 = a1 + 2sh;
a4 = -a1;

pair uv[];
numeric len;
len = .8gl;
uv0 = .38[a1,a2] + .8gw*down;
uv1 = uv0 + len*right;
uv2 = uv0 + len*up;

%% reference frame
drawoptions(withcolor DeepRed);
drawarrow uv0 -- uv1;
drawarrow uv0 -- uv2;
drawoptions();

%% arrows
interim ahlength := 7bp;
interim ahangle := 20;

drawarrow a1 -- (a1 + .5sh);
drawarrow a2 -- (a2 + .5sh);
drawarrow a3{down} ... (.7[a1,a2] + 1.3gw*down) ... (a4+1down){up};

%% labels
text_directive := "\footnotesize";
dotlabel.lft(textext(mytext("g")), a1);
dotlabel.lft(textext(mytext("g")), a2);
dotlabel.lft(textext(mytext("g")), a3);

endfig;


end;
