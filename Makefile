graphics = \
  illustrations/charts-for-element.eps          illustrations/misalign-fig4.eps\
  illustrations/connecting-blocks.eps           illustrations/misalign-fig5.eps\
  illustrations/DNA-database-concept.eps        illustrations/misalign-fig6.eps\
  illustrations/DNA-database-tutorial.eps       illustrations/misaligning-element.eps\
  illustrations/fibre.eps                       illustrations/misalignment-for-element.eps\
  illustrations/formula1.eps                    illustrations/model-cell.eps\
  illustrations/formula2.eps                    illustrations/model-fig8.eps\
  illustrations/formula4.eps                    illustrations/model-JLAB-DNA.eps\
  illustrations/formula5.eps                    illustrations/model-JLAB.eps\
  illustrations/geo-routines-1.eps              illustrations/model-JLAB-fibres.eps\
  illustrations/geo-routines-2.eps              illustrations/model-PSR.eps\
  illustrations/integration-node-and-fibre.eps  illustrations/model-RHIC.eps\
  illustrations/integration-nodes.eps           illustrations/model-simple.eps\
  illustrations/LEGO-concept.eps                illustrations/model-tutorial.eps\
  illustrations/LEGO-element-ref-frame.eps      illustrations/patching-beam-lines.eps\
  illustrations/LEGO.eps                        illustrations/patching-element.eps\
  illustrations/LEGO-faces.eps                  illustrations/pseudo-Euclidean-maps.eps\
  illustrations/misaligned-planar-fibre.eps     illustrations/rbend.eps\
  illustrations/misalign-fig1.eps               illustrations/recutDKD.eps\
  illustrations/misalign-fig2.eps               illustrations/recutMKM.eps\
  illustrations/misalign-fig3.eps               illustrations/space-charge-kick.eps

texfiles = PTC-LibUG.tex \
  front/titles.tex \
  front/copyright.tex \
  chapters/chap01.tex \
  chapters/chap02.tex \
  chapters/chap03.tex \
  chapters/chap04.tex \
  chapters/chap05.tex \
  chapters/chap06.tex \
  chapters/chap07.tex \
  chapters/chap08.tex \
  chapters/chap09.tex \
  chapters/chap10.tex \
  chapters/chap11.tex \
  chapters/gloss.tex \
  appends/appenA.tex \
  appends/appenB.tex \
  appends/appenC.tex \
  appends/appenD.tex \
  appends/appenE.tex \
  appends/appenF.tex

cfgfiles = \
  PTC-LibUG.ist \
  PTC-LibUG.gst \
  ptccmds.gst

bibfiles = PTCbibliography.bib

default: $(texfiles) $(cfgfiles) $(bibfiles)
	pdflatex PTC-LibUG
	bibtex PTC-LibUG
	makeindex -s PTC-LibUG.ist -t PTC-LibUG.ilg -o PTC-LibUG.ind PTC-LibUG.idx
	makeindex -s PTC-LibUG.gst -t PTC-LibUG.glg -o PTC-LibUG.gls PTC-LibUG.glo
	makeindex -s ptccmds.gst -t ptccmds.glg -o ptccmds.gls ptccmds.glo
	pdflatex PTC-LibUG
	pdflatex PTC-LibUG

clean:
	rm -f *~ *.blg *.glg *.ilg *.log *.out
	rm -f front/*~ chapters/*~ appends/*~

cleanall:
	rm -f *~ *.aux *.bbl *.blg *.glg *.glo *.gls *.idx *.ilg *.ind *.lof *.log *.lot *.out *.pdf *.toc
	rm -f front/*~ front/*.aux
	rm -f chapters/*~ chapters/*.aux
	rm -f appends/*~ appends/*.aux

