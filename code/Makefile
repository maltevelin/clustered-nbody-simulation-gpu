FUTHARKOPENCL=futhark opencl
FUTHARKC=futhark c --Werror --safe
FUTHARKBENCH=futhark bench --backend=opencl
FUTHARKAUTOTUNE=futhark autotune --backend=opencl

TESTBODIES=$(patsubst %,data/%-bodies.in,100 1000 10000 100000 200000 1000000 2000000 4000000)

BENCHES=nbody n_squared_body

.PHONY: run clean benchmarks test

.PRECIOUS: data/%-bodies.in

all: nbody.py

test:
	futhark test --exclude=notest --backend=opencl --compiled nbody.fut # n_squared_body.fut


nbody.py: nbody.fut
	futhark pyopencl --library nbody.fut

nbody.c: nbody.fut
	futhark opencl --library nbody.fut

_nbody.so: nbody.c
	build_futhark_ffi nbody

nbody: nbody.fut
	$(FUTHARKOPENCL) $<

n_squared_body: n_squared_body.fut
	$(FUTHARKOPENCL) $<

benchmarks: $(patsubst %,bench_results/%.txt,$(BENCHES))

%.tuning: %.fut
	$(FUTHARKAUTOTUNE) $<

data/%-bodies.in: data/gen_bin_input.sh
	@mkdir -p $(@D)
	$(@D)/gen_bin_input.sh $* 1 0.5 > $@

bench_results/%.txt: %.fut $(TESTBODIES)
	@mkdir -p $(@D)
	$(FUTHARKBENCH) --no-tuning $< | tee $@

bench_results/%_tuned.txt: %.fut $(TESTBODIES) %.tuning
	@mkdir -p $(@D)
	$(FUTHARKBENCH) $< | tee $@

clean:
	rm -rf nbody.py nbody nnbody *.pyc *.so *.o *.c *.tuning __pycache__ bench_results

run: nbody.py
	@python3 nbody-gui.py
