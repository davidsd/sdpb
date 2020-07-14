CONTAINERS=sdpb.sif base.sif
SINGULARITY=singularity
SINGULARITY_ROOT=sudo $(SINGULARITY)

.PHONY: all
all: sdpb.sif sdpb_tested.txt

.PHONY: clean clean_tests
clean_tests:
	rm -rf test
	rm sdpb_tested.txt

clean: clean_tests
	rm -f ${CONTAINERS}

%.sif: %.def
	$(SINGULARITY_ROOT) build --force --notest $@ $<

sdpb.sif: base.sif $(shell find ../src -type f)

sdpb_tested.txt: sdpb.sif
	$(SINGULARITY) test sdpb.sif&&touch sdpb_tested.txt
