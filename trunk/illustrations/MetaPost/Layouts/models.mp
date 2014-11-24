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


%:figure 0: basic cells
beginfig(0)
numeric uv;
uv := .3cm;

picture cell[];
numeric scl;
scl := 5;

%% cell with rectangular bend
drawoptions(withpen pencircle scaled dwt);
draw ((rm * dir(90 - .5da)) -- (rp * dir(90)) -- (rm * dir(90 + .5da)))
  scaled scl;
drawoptions(withpen pencircle scaled pw withcolor Crimson);
fill (rect(lb,ew) shifted (rho*dir(90))) scaled scl withcolor gray;
draw (rect(lb,ew) shifted (rho*dir(90))) scaled scl;
draw (subpath (t1,t2) of circ0) shifted (rho*dir(90)) scaled scl dashed pdots;
drawoptions(withpen pencircle scaled pw withcolor DarkBlue);
fill (rect(lq,ew) rotated ( .5da) shifted (rq*dir(90 + aq))) scaled scl
  withcolor gray;
draw (rect(lq,ew) rotated ( .5da) shifted (rq*dir(90 + aq))) scaled scl;
fill (rect(lq,ew) rotated (-.5da) shifted (rq*dir(90 - aq))) scaled scl
  withcolor gray;
draw (rect(lq,ew) rotated (-.5da) shifted (rq*dir(90 - aq))) scaled scl;
drawoptions(withpen pencircle scaled 2pw);
drawuvecs(uv, scl*(.5(-lb,-sag) shifted (rho*dir(90))),  .5da, Crimson);
drawuvecs(uv, scl*(.5( lb,-sag) shifted (rho*dir(90))), -.5da, Crimson);
drawuvecs(uv, scl*(rq*dir(90 + aq) shifted ( .5lq*dir(180 + .5da))),  .5da, DarkBlue);
drawuvecs(uv, scl*(rq*dir(90 + aq) shifted (-.5lq*dir(180 + .5da))),  .5da, DarkBlue);
drawuvecs(uv, scl*(rq*dir(90 - aq) shifted ( .5lq*dir(180 - .5da))), -.5da, DarkBlue);
drawuvecs(uv, scl*(rq*dir(90 - aq) shifted (-.5lq*dir(180 - .5da))), -.5da, DarkBlue);
drawoptions();
cell1 = currentpicture;
currentpicture := nullpicture;

%% cell with ``straight'' bend
drawoptions(withpen pencircle scaled dwt);
draw ((rm * dir(90 - .5da)) -- (rp * dir(90)) -- (rm * dir(90 + .5da)))
  scaled scl;
drawoptions(withpen pencircle scaled pw withcolor Hunter);
fill (rect(lb,ew) shifted ((rho-.5sag)*dir(90))) scaled scl withcolor gray;
draw (rect(lb,ew) shifted ((rho-.5sag)*dir(90))) scaled scl;
draw ((-.5lb,-.5sag) -- (.5lb,-.5sag)) shifted (rho*dir(90)) scaled scl dashed pdots;
drawoptions(withpen pencircle scaled pw withcolor DarkBlue);
fill (rect(lq,ew) rotated ( .5da) shifted (rq*dir(90 + aq))) scaled scl
  withcolor gray;
draw (rect(lq,ew) rotated ( .5da) shifted (rq*dir(90 + aq))) scaled scl;
fill (rect(lq,ew) rotated (-.5da) shifted (rq*dir(90 - aq))) scaled scl
  withcolor gray;
draw (rect(lq,ew) rotated (-.5da) shifted (rq*dir(90 - aq))) scaled scl;
drawoptions(withpen pencircle scaled 2pw);
drawuvecs(uv, scl*(.5(-lb,0) shifted ((rho-.5sag)*dir(90))), 0, Hunter);
drawuvecs(uv, scl*(.5( lb,0) shifted ((rho-.5sag)*dir(90))), 0, Hunter);
drawuvecs(uv, scl*(rq*dir(90 + aq) shifted ( .5lq*dir(180 + .5da))),  .5da, DarkBlue);
drawuvecs(uv, scl*(rq*dir(90 + aq) shifted (-.5lq*dir(180 + .5da))),  .5da, DarkBlue);
drawuvecs(uv, scl*(rq*dir(90 - aq) shifted ( .5lq*dir(180 - .5da))), -.5da, DarkBlue);
drawuvecs(uv, scl*(rq*dir(90 - aq) shifted (-.5lq*dir(180 - .5da))), -.5da, DarkBlue);
drawoptions();
cell2 = currentpicture;
currentpicture := nullpicture;

draw cell1;
draw cell2 shifted (.72rm * scl * right);

%% labels
text_directive := "\footnotesize";
label.bot(textext(mytext("cell for figure eight and collider")),
  scl*(.9rho)*dir(90));
label.bot(textext(mytext("cell for ring with forward and reverse propagation")),
  scl*(.9rho)*dir(90)) shifted (.72rm * scl * right);

endfig;


%:figure 1: accel building blocks (DNA layouts)
beginfig(1)

%% global frame
drawoptions(withcolor Hunter);
drawuvecs(1.2u, p5, 0, Hunter);
fill fullcircle scaled 2.5 shifted p5;
draw fullcircle scaled 4.0 shifted p5;
drawoptions();

%% figure eight
draw psr shifted p1;
draw psr_minus shifted (p1 + p4);
drawoptions(withcolor Crimson dashed pdots);
draw arcp rotated -.5da shifted (p1 + rho*dir(2da));
draw arcp rotated  .5da shifted (p1 + rho*dir(3da));
drawoptions();
draw arrcw rotated 45 shifted p1;
draw arrcw rotated -135 shifted p1;
draw arrccw rotated 45 shifted (p1 + p4);
draw arrccw rotated -135 shifted (p1 + p4);
%% concentric PSRs
draw psr;
draw psr_qfb;
draw arrcw;
draw arrccw rotated 180;
%% collider
draw psr shifted p2;
draw psr_minus shifted (p2 + p4);
draw arcp rotated -.5da shifted (p2 + rho*dir(2da));
draw arcp rotated  .5da shifted (p2 + rho*dir(3da));
draw arrcw rotated 45 shifted p2;
draw arrcw rotated -135 shifted p2;
draw arrcw rotated 45 shifted (p2 + p4);
draw arrcw rotated -135 shifted (p2 + p4);

%% labels
text_directive := "\footnotesize";
label(textext(mytext("figure eight (\texttt{Fig8})")), p3 + p1);
label(textext(mytext("forward and reverse rings (\texttt{PSR1} and \texttt{PSR2})")), p3);
label(textext(mytext("collider (\texttt{Col1} and \texttt{Col2})")), p3 + p2);
drawoptions(withcolor DarkCyan);
label(textext(mytext("\texttt{Fig8}")), p1);
label(textext(mytext("\texttt{PSR1}")),  .6p5);
label(textext(mytext("\texttt{PSR2}")), -.6p5);
label(textext(mytext("\texttt{Col1}")), p2);
label(textext(mytext("\texttt{Col2}")), p2 + p4);
drawoptions();
drawoptions(withcolor Hunter);
label.rt(textext(mytext(
  "\begin{minipage}{.5in} "
  & "global ref.\\ frame"
  & "\end{minipage}")),
  origin) shifted (p5 + 18pt * up) withcolor Hunter;
label.top(textext(mytext("$X$")), p5 + 1.2u * up);
label.bot(textext(mytext("$Y$")), p5);
label.rt(textext(mytext("$Z$")), p5 + 1.2u * right);
drawoptions();

endfig;


%:figure 2: accel building blocks (DNA layouts)
beginfig(2)
draw psr shifted p1;
draw psr_minus rotated 180;
draw psr_qfb shifted p2;

text_directive := "\footnotesize\texttt";
label(textext(mytext("\texttt{build\_PSR}")), p3 + p1);
label(textext(mytext("\texttt{build\_PSR\_minus}")), p3);
label(textext(mytext("\texttt{build\_Quad\_for\_Bend}")), p3 + p2);

endfig;


%:figure 3: fig8 w/ text re fibres
beginfig(3)

draw psr;
draw psr_minus shifted p4;
drawoptions(withcolor Crimson dashed pdots);
draw arcp rotated -.5da shifted (rho*dir(2da));
draw arcp rotated  .5da shifted (rho*dir(3da));
drawoptions();
draw arrcw rotated 45;
draw arrcw rotated -135;
draw arrccw rotated 45 shifted p4;
draw arrccw rotated -135 shifted p4;

%% highlight common fibres
(p7 - p5) dotprod up = 0;
(p7 - (rho - .5sag) * dir(90 - .5da)) dotprod dir(90 - .5da) = 0;
p8 = p7 reflectedabout(origin, p5);
linecap := butt;
draw p7 -- p8 withpen pencircle scaled dwt withcolor Green;
linecap := rounded;
drawoptions(withpen pencircle scaled pw withcolor Green);
fill rect(lq,ew) shifted (rq*dir(90 -.5da + aq)) withcolor gray;
draw rect(lq,ew) shifted (rq*dir(90 -.5da + aq));
fill rect(lq,ew) shifted (rq*dir(90 +.5da - aq)) withcolor gray;
draw rect(lq,ew) shifted (rq*dir(90 +.5da - aq));
drawoptions();

%% labels
text_directive := "\footnotesize";
drawoptions(withcolor DarkCyan);
label(textext(mytext("\texttt{Fig8}")), origin);
drawoptions();
text_directive := "\small";
label.rt(textext(mytext(%
  "\begin{minipage}{4.5cm}\RaggedRight\noindent"
  & "64 fibres in layout \texttt{fig8} point to "
  & "64 fibres in DNA layout \texttt{L3}"
  & "\end{minipage}"
)), .5p2 + p4);
label.rt(textext(mytext(%
  "\begin{minipage}{4.5cm}\RaggedRight\noindent"
  & "12 fibres in layout \texttt{fig8} point to "
  & "6 fibres in DNA layout \texttt{L4}"
  & "\end{minipage}"
)), .5p2 + p5) withcolor Green;
label.rt(textext(mytext(%
  "\begin{minipage}{4.5cm}\RaggedRight\noindent"
  & "64 fibres in layout \texttt{fig8} point to "
  & "64 fibres in DNA layout \texttt{L4}"
  & "\end{minipage}"
)), .5p2);

endfig;


%:figure 4: geometry of rectangular bend
beginfig(4)
numeric uv;
uv := .5cm;

picture cell[];
numeric scl;
scl := 6;

%% cell with rectangular bend
drawoptions(withpen pencircle scaled pw withcolor Crimson);
fill rect(lb,ew) scaled scl withcolor gray;
draw rect(lb,ew) scaled scl;
draw (subpath (t1,t2) of circ0) scaled scl dashed pdots;

drawoptions(withpen pencircle scaled pw);
draw (.5(-lb,-sag) -- ctr1) scaled scl dashed evenly;
draw (.5( lb,-sag) -- ctr1) scaled scl dashed evenly;
drawoptions(withpen pencircle scaled .5pw);
draw (.5(-lb,-sag) -- .5( lb,-sag)) scaled scl;

drawoptions(withpen pencircle scaled 2pw);
drawuvecs(uv, scl * .5(-lb,-sag),  .5da, Crimson);
drawuvecs(uv, scl * .5( lb,-sag), -.5da, Crimson);
drawoptions();

drawarrow (.3[scl*ctr1, .5scl*(-lb,-sag)]{dir( .5da)}
        .. .3[scl*ctr1, .5scl*( lb,-sag)]{dir(-.5da)});

%% labels
text_directive := "\footnotesize\texttt";
label.top(textext(mytext("ang")), .3[scl*ctr1,origin]);
label.top(textext(mytext("Larc = 2.54948")), 3*up);

endfig;


%:figure 5: small PSR1  and PSR2
beginfig(5)

%% concentric PSRs
draw psr scaled .85;
draw psr_qfb  scaled .85;
draw arrcw scaled .85;
draw arrccw rotated 180 scaled .85;

%% labels
text_directive := "\footnotesize";
drawoptions(withcolor DarkCyan);
label(textext(mytext("\texttt{PSR1}")),  .6p5 scaled .80);
label(textext(mytext("\texttt{PSR2}")), -.6p5 scaled .80);
drawoptions();

endfig;


%:figure 6: matching rings
beginfig(6)

%% rings
draw psr shifted p1 scaled .85;
draw psr_minus shifted (p1 + p4 + 25up) scaled .85;

%% arrow
interim ahlength := 5bp;
interim ahangle := 15;
drawarrow ((p1 + (rm + 25)*up + 11*right) -- (p1 + rm*up + 11*right)) scaled 0.85 dashed pdots withcolor .4white;

%% labels
text_directive := "\footnotesize";
drawoptions(withcolor DarkCyan);
label(textext(mytext("\texttt{L4}")), p1 scaled .85);
label(textext(mytext("\texttt{L3}")), (p1 + p4 + 25up) scaled .85);
drawoptions();

endfig;


%:figure 7: small Col1 and Col2
beginfig(7)

%% collider
draw psr shifted p2 scaled .85;
draw psr_minus shifted (p2 + p4) scaled .85;
draw arcp rotated -.5da shifted (p2 + rho*dir(2da)) scaled .85;
draw arcp rotated  .5da shifted (p2 + rho*dir(3da)) scaled .85;
draw arrcw rotated 45 shifted p2 scaled .85;
draw arrcw rotated -135 shifted p2 scaled .85;
draw arrcw rotated 45 shifted (p2 + p4) scaled .85;
draw arrcw rotated -135 shifted (p2 + p4) scaled .85;

%% labels
text_directive := "\footnotesize";
drawoptions(withcolor DarkCyan);
label(textext(mytext("\texttt{Col1}")), p2 scaled .85);
label(textext(mytext("\texttt{Col2}")), (p2 + p4) scaled .85);
drawoptions();

endfig;


%:figure 8: zoom of Col1 and Col2
beginfig(8)
numeric myscale;
myscale := 1.22;

%% collider
draw psr shifted p2 scaled myscale;
draw psr_minus shifted (p2 + p4) scaled myscale;
draw arcp rotated -.5da shifted (p2 + rho*dir(2da)) scaled myscale;
draw arcp rotated  .5da shifted (p2 + rho*dir(3da)) scaled myscale;
draw arrcw rotated 45 shifted p2 scaled myscale;
draw arrcw rotated -135 shifted p2 scaled myscale;
drawoptions(withpen pencircle scaled pw withcolor FireBrick);
draw arrcw rotated 45 shifted (p2 + p4) scaled myscale;
draw arrcw rotated -135 shifted (p2 + p4) scaled myscale;
drawoptions();

%% clipping
path myzoom;
myzoom := rect(1.50rp,1.13rp) shifted (.92p2 + .5p4) scaled myscale;
clip currentpicture to myzoom;
%draw myzoom;

%% labels
numeric ld;
ld := arclength(rm*up -- rq*dir(90 - .5da + aq));

text_directive := "\footnotesize";
drawoptions(withcolor DarkCyan);
label(textext(mytext("\texttt{Col1}")), (p2 + .3p4) scaled myscale);
drawoptions(withcolor FireBrick);
label(textext(mytext("\texttt{Col2}")), (p2 + .7p4) scaled myscale);

%text_directive := "\scriptsize";
text_directive := "\footnotesize";
drawoptions(withcolor DarkCyan);
label(textext(mytext("1")),  (p2 + .88rm*up + .5*(ld - .5lq)*right) scaled myscale);
label(textext(mytext("70")), (p2 + .88rm*up + .5*(ld - .5lq)*left)  scaled myscale);
label(textext(mytext("4")),  (p2 + .82rp*dir(90 - .5da)) scaled myscale);
label(textext(mytext("67")), (p2 + .82rp*dir(90 + .5da)) scaled myscale);
drawoptions(withcolor FireBrick);
label(textext(mytext("3")),  (p2+p4 + .88rm*down + .5*(ld - .5lq)*right) scaled myscale);
label(textext(mytext("4")),  (p2+p4 + .88rm*down + .5*(ld - .5lq)*left)  scaled myscale);
label(textext(mytext("7")),  (p2+p4 + .82rp*dir(-90 - .5da))  scaled myscale);
label(textext(mytext("70")), (p2+p4 + .82rp*dir(-90 + .5da))  scaled myscale);
label(textext(mytext("14")), (p2+p4 + .82rp*dir(-90 - 1.5da))  scaled myscale);
drawoptions();

endfig;


%:figure 9: second zoom of Col1 and Col2
beginfig(9)
numeric myscale;
myscale := 1.22;

%% collider
draw psr shifted p2 scaled myscale;
draw psr_minus shifted (p2 + p4) scaled myscale;
draw arcp rotated -.5da shifted (p2 + rho*dir(2da)) scaled myscale;
draw arcp rotated  .5da shifted (p2 + rho*dir(3da)) scaled myscale;
draw arrcw rotated 45 shifted p2 scaled myscale;
draw arrcw rotated -135 shifted p2 scaled myscale;
drawoptions(withpen pencircle scaled pw withcolor FireBrick);
draw arrcw rotated 45 shifted (p2 + p4) scaled myscale;
draw arrcw rotated -135 shifted (p2 + p4) scaled myscale;
drawoptions();

%% clipping
path myzoom;
myzoom := rect(1.2rp,0.9rp) shifted (p2 + .5p4) scaled myscale;
clip currentpicture to myzoom;
%draw myzoom;

%% cross-hairs
numeric chlen; chlen := 2.5pt;
picture crosshair;
crosshair := image(
  drawoptions(withpen pencircle scaled .5pw withcolor black);
  draw (-chlen, 0) -- (chlen, 0);
  draw (0, -chlen) -- (0, chlen);
  drawoptions();
);
draw crosshair shifted ((p2 + .5p4 + 1.7lb*right) scaled myscale);
draw crosshair shifted ((p2 + .5p4 - 1.7lb*right) scaled myscale);

%% labels
text_directive := "\footnotesize";
drawoptions(withcolor DarkCyan);
label(textext(mytext("\texttt{Col1}")), (p2 + .3p4) scaled myscale);
drawoptions(withcolor FireBrick);
label(textext(mytext("\texttt{Col2}")), (p2 + .7p4) scaled myscale);

endfig;


end;
