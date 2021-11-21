.PHONY: clean
clean:
	rm -f *.txt

.PHONY: run
run:
	julia 2-julia_sho.jl

.PHONY: total
total: clean run
