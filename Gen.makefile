NAB=$(AMBERHOME)/bin/nab
LEAP=$(AMBERHOME)/bin/tleap

RNAPDBS = gcaucg.pdb gcgccg.pdb gcgucg.pdb gugucg.pdb gcggcg.pdb gcaacg.pdb guggug.pdb

all: $(subst .pdb,_GMX.top,$(RNAPDBS))

%.nabbin: %.nab
	$(NAB) -o $@ $<

$(RNAPDBS): genrnaah.nabbin
	./$<

%.leap: base.leap.in
	sed -e "s/%S%/$*/g" $< > $@

%.ambtop %.crd: %.leap %.pdb
	$(LEAP) -f $<

%_GMX.top %_GMX.gro: %.ambtop %.crd
	./acpype.py -p $*.ambtop -x $*.crd

acpype.py: ../acpype/acpype.py
	ln -s $< $@


