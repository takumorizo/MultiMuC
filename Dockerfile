FROM julia:1.2

RUN set -x; \
	apt-get update; \
	apt-get install -y sudo; \
	sudo apt-get install -y git; \
	cd /opt; \
	git clone https://github.com/takumorizo/MultiMuC.git; \
	cd /opt/MultiMuC/src; \
	julia --project -e "using Pkg; Pkg.instantiate();"; \
	mkdir -p /opt/MultiMuC/test; \
	julia --project ./main.jl MAP --seed 0 -n 600 -b 100 \
	--pdata=0.02 --pevidence=0.2 --pcon=0.999 --theta=0.5 --rho=0.1 \
	-o /opt/MultiMuC/test /opt/MultiMuC/examples/input/Strelka2/Strelka2_HCC1954_tree_3_med_L_a_0.txt	/opt/MultiMuC/examples/input/Strelka2/Strelka2_HCC1954_tree_3_med_H.txt; \
	cd /opt/MultiMuC/test; \
	md5sum *.txt > ./check; \
	cat ./check; \
	cd /opt/MultiMuC/examples/output/Strelka2/; \
	md5sum -c /opt/MultiMuC/test/check;

ENTRYPOINT ["julia --project /opt/MultiMuC/src/main.jl"]
